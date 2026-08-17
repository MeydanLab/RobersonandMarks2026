#!/usr/bin/env python3
"""
Quantify C->G mismatches in paired-end RNA-seq BAMs by transcript feature and CDS codon phase.

Example:
    python3 rnaseq_c2g_feature_codon_analysis.py \
        --bam-dir /path/to/bam_directory \
        --bam-glob '*_Aligned.sortedByCoord.out.bam' \
        --gtf /path/to/annotation.gtf \
        --fasta /path/to/reference.fa \
        --output-prefix dataset_C2G_feature_codon \
        --samtools samtools \
        --samtools-threads 4

Key rules implemented here:
1. BAMs are streamed through `samtools view`; no third-party Python packages are required.
2. Only bases actually sequenced/aligned in each read are counted. Spliced-out introns (`N` in CIGAR)
   do not contribute any coverage because the read does not cover those bases.
3. Both mates are analyzed independently as aligned reads. No fragment gap is inferred between mates.
4. Feature and codon-phase assignments come from the GTF. If overlapping coding transcripts disagree
   on a covered genomic position, that covered reference C is marked ambiguous and excluded from the
   feature/codon denominators.
5. Rates are computed exactly as requested: C->G mismatches / covered reference C bases.
6. Counting is read-based, not locus-collapsed: if multiple reads cover the same genomic C, each observed
   C contributes to the denominator, and each observed C->G event contributes to the numerator.

Outputs:
- <prefix>.dataset_summary.tsv
- <prefix>.feature_summary.tsv
- <prefix>.codon_summary.tsv
- <prefix>.run_metadata.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import pickle
import posixpath
import shutil
import subprocess
import sys
import zipfile
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Sequence, Tuple
from xml.etree import ElementTree as ET


BASES = {"A", "C", "G", "T"}
CIGAR_OPS = set("MIDNSHP=X")
DEFAULT_SKIP_FLAG = 0x4 | 0x100 | 0x800
BUCKET_SIZE = 100_000
FEATURE_ORDER = ["5' UTR", "CDS", "3' UTR", "intron"]
CODON_ORDER = [1, 2, 3]
CACHE_VERSION = 1
ANNOT_UNANNOTATED = 0
ANNOT_AMBIGUOUS = 1
ANNOT_5UTR = 2
ANNOT_CDS_1 = 3
ANNOT_CDS_2 = 4
ANNOT_CDS_3 = 5
ANNOT_CDS_AMBIGUOUS = 6
ANNOT_3UTR = 7
ANNOT_INTRON = 8
XLSX_MAIN_NS = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
XLSX_REL_NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
XLSX_PKG_REL_NS = {"pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships"}


@dataclass(frozen=True)
class FastaIndexEntry:
    length: int
    offset: int
    bases_per_line: int
    line_width: int


class IndexedFasta:
    def __init__(self, fasta_path: Path):
        self.fasta_path = fasta_path
        self.fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
        if not self.fai_path.exists():
            raise FileNotFoundError(
                f"Missing FASTA index: {self.fai_path}. Create it first with `samtools faidx`."
            )
        self.index = self._load_index(self.fai_path)
        self._fh = fasta_path.open("rb")
        self._cache_name = ""
        self._cache_start = 0
        self._cache_end = 0
        self._cache_seq = ""

    @staticmethod
    def _load_index(path: Path) -> Dict[str, FastaIndexEntry]:
        index: Dict[str, FastaIndexEntry] = {}
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                index[fields[0]] = FastaIndexEntry(
                    length=int(fields[1]),
                    offset=int(fields[2]),
                    bases_per_line=int(fields[3]),
                    line_width=int(fields[4]),
                )
        return index

    def close(self) -> None:
        self._fh.close()

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if start0 < 0 or end0 < start0:
            raise ValueError(f"Invalid FASTA range: {chrom}:{start0}-{end0}")
        if start0 == end0:
            return ""
        if (
            chrom == self._cache_name
            and start0 >= self._cache_start
            and end0 <= self._cache_end
        ):
            rel_start = start0 - self._cache_start
            rel_end = end0 - self._cache_start
            return self._cache_seq[rel_start:rel_end]

        entry = self.index.get(chrom)
        if entry is None:
            raise KeyError(f"Reference sequence `{chrom}` not found in {self.fasta_path}")
        if end0 > entry.length:
            raise ValueError(f"Requested {chrom}:{start0}-{end0} past contig length {entry.length}")

        seq = self._fetch_uncached(start0, end0, entry)
        span = end0 - start0
        if span < 1_000_000:
            cache_end = min(entry.length, max(end0, start0 + 1_000_000))
            self._cache_name = chrom
            self._cache_start = start0
            self._cache_end = cache_end
            self._cache_seq = self._fetch_uncached(start0, cache_end, entry)
        return seq

    def _fetch_uncached(self, start0: int, end0: int, entry: FastaIndexEntry) -> str:
        start_line = start0 // entry.bases_per_line
        start_col = start0 % entry.bases_per_line
        end_line = (end0 - 1) // entry.bases_per_line
        end_col = (end0 - 1) % entry.bases_per_line

        byte_start = entry.offset + start_line * entry.line_width + start_col
        byte_end = entry.offset + end_line * entry.line_width + end_col + 1

        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()


def parse_cigar(cigar: str) -> Iterator[Tuple[int, str]]:
    digits: List[str] = []
    for char in cigar:
        if char.isdigit():
            digits.append(char)
            continue
        if char not in CIGAR_OPS or not digits:
            raise ValueError(f"Invalid CIGAR: {cigar}")
        yield int("".join(digits)), char
        digits.clear()
    if digits:
        raise ValueError(f"Trailing length in CIGAR: {cigar}")


def aligned_blocks(cigar: str, ref_start0: int) -> Iterator[Tuple[int, int, int]]:
    """
    Yield (read_start0, ref_start0, length) blocks for aligned query/reference sequence.
    """
    read_pos = 0
    ref_pos = ref_start0
    for length, op in parse_cigar(cigar):
        if op in {"M", "=", "X"}:
            yield read_pos, ref_pos, length
            read_pos += length
            ref_pos += length
        elif op in {"I", "S"}:
            read_pos += length
        elif op in {"D", "N"}:
            ref_pos += length
        elif op in {"H", "P"}:
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation `{op}` in {cigar}")


def parse_gtf_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        if not key or not value:
            continue
        attrs[key] = value.strip().strip('"')
    return attrs


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return []
    merged = [sorted_intervals[0]]
    for start0, end0 in sorted_intervals[1:]:
        last_start0, last_end0 = merged[-1]
        if start0 <= last_end0:
            merged[-1] = (last_start0, max(last_end0, end0))
        else:
            merged.append((start0, end0))
    return merged


def subtract_intervals(
    source: Sequence[Tuple[int, int]],
    subtractors: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    remaining: List[Tuple[int, int]] = []
    sub_idx = 0
    subtractors = sorted(subtractors)
    for src_start0, src_end0 in source:
        cursor = src_start0
        while sub_idx < len(subtractors) and subtractors[sub_idx][1] <= src_start0:
            sub_idx += 1
        scan_idx = sub_idx
        while scan_idx < len(subtractors):
            sub_start0, sub_end0 = subtractors[scan_idx]
            if sub_start0 >= src_end0:
                break
            if sub_start0 > cursor:
                remaining.append((cursor, min(sub_start0, src_end0)))
            cursor = max(cursor, sub_end0)
            if cursor >= src_end0:
                break
            scan_idx += 1
        if cursor < src_end0:
            remaining.append((cursor, src_end0))
    return [interval for interval in remaining if interval[0] < interval[1]]


def interval_contains(
    intervals: Sequence[Tuple[int, int]],
    starts: Sequence[int],
    pos0: int,
) -> int | None:
    idx = bisect_right(starts, pos0) - 1
    if idx < 0:
        return None
    start0, end0 = intervals[idx]
    if start0 <= pos0 < end0:
        return idx
    return None


@dataclass(frozen=True)
class CDSChunk:
    start0: int
    end0: int
    cds_offset0: int


@dataclass
class TranscriptRecord:
    transcript_id: str
    gene_id: str = ""
    gene_name: str = ""
    chrom: str = ""
    strand: str = ""
    source: str = ""
    tx_start0: int | None = None
    tx_end0: int | None = None
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    source: str
    tx_start0: int
    tx_end0: int
    exons: List[Tuple[int, int]]
    exon_starts: List[int]
    introns: List[Tuple[int, int]]
    intron_starts: List[int]
    cds_intervals: List[Tuple[int, int]]
    cds_starts: List[int]
    cds_chunks: List[CDSChunk]
    cds_min0: int
    cds_max0: int

    def annotate_position(self, pos0: int) -> Tuple[str | None, int | None]:
        cds_idx = interval_contains(self.cds_intervals, self.cds_starts, pos0)
        if cds_idx is not None:
            chunk = self.cds_chunks[cds_idx]
            if self.strand == "+":
                cds_offset0 = chunk.cds_offset0 + (pos0 - chunk.start0)
            else:
                cds_offset0 = chunk.cds_offset0 + (chunk.end0 - 1 - pos0)
            return "CDS", (cds_offset0 % 3) + 1

        exon_idx = interval_contains(self.exons, self.exon_starts, pos0)
        if exon_idx is not None:
            if self.strand == "+":
                feature = "5' UTR" if pos0 < self.cds_min0 else "3' UTR"
            else:
                feature = "5' UTR" if pos0 >= self.cds_max0 else "3' UTR"
            return feature, None

        intron_idx = interval_contains(self.introns, self.intron_starts, pos0)
        if intron_idx is not None:
            return "intron", None

        return None, None


@dataclass
class DatasetStats:
    reads_seen: int = 0
    reads_skipped_nm_zero: int = 0
    total_c_bases: int = 0
    total_c_to_g: int = 0
    feature_c_bases: Counter = field(default_factory=Counter)
    feature_c_to_g: Counter = field(default_factory=Counter)
    codon_c_bases: Counter = field(default_factory=Counter)
    codon_c_to_g: Counter = field(default_factory=Counter)
    ambiguous_feature_c_bases: int = 0
    ambiguous_feature_c_to_g: int = 0
    unannotated_c_bases: int = 0
    unannotated_c_to_g: int = 0
    ambiguous_codon_c_bases: int = 0
    ambiguous_codon_c_to_g: int = 0
    missing_reference_reads: int = 0

    def add(self, other: "DatasetStats") -> None:
        self.reads_seen += other.reads_seen
        self.reads_skipped_nm_zero += other.reads_skipped_nm_zero
        self.total_c_bases += other.total_c_bases
        self.total_c_to_g += other.total_c_to_g
        self.feature_c_bases.update(other.feature_c_bases)
        self.feature_c_to_g.update(other.feature_c_to_g)
        self.codon_c_bases.update(other.codon_c_bases)
        self.codon_c_to_g.update(other.codon_c_to_g)
        self.ambiguous_feature_c_bases += other.ambiguous_feature_c_bases
        self.ambiguous_feature_c_to_g += other.ambiguous_feature_c_to_g
        self.unannotated_c_bases += other.unannotated_c_bases
        self.unannotated_c_to_g += other.unannotated_c_to_g
        self.ambiguous_codon_c_bases += other.ambiguous_codon_c_bases
        self.ambiguous_codon_c_to_g += other.ambiguous_codon_c_to_g
        self.missing_reference_reads += other.missing_reference_reads


@dataclass
class CPositionBucket:
    positions_by_code: List[List[int]]
    code_by_position: Dict[int, int]


def normalize_label(value: str | None) -> str:
    return value if value else ""


def normalize_gene_key(value: str | None) -> str:
    return normalize_label(value).strip()


def excel_column_index(cell_ref: str) -> int:
    column_text: List[str] = []
    for char in cell_ref:
        if char.isalpha():
            column_text.append(char.upper())
        else:
            break
    index = 0
    for char in column_text:
        index = (index * 26) + (ord(char) - ord("A") + 1)
    return max(index - 1, 0)


def load_allowed_genes_from_workbook(
    workbook_path: Path,
    dataset_label: str | None,
    minimum_ratio: float,
) -> set[str]:
    with zipfile.ZipFile(workbook_path) as archive:
        shared_strings: List[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall("main:si", XLSX_MAIN_NS):
                shared_strings.append("".join(item.itertext()))

        workbook_root = ET.fromstring(archive.read("xl/workbook.xml"))
        rel_root = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rel_root.findall("pkgrel:Relationship", XLSX_PKG_REL_NS)
            if "Id" in rel.attrib and "Target" in rel.attrib
        }

        target_path = None
        for sheet in workbook_root.findall("main:sheets/main:sheet", XLSX_REL_NS):
            if sheet.attrib.get("name") != "All_gene_pass_rates":
                continue
            rel_id = sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
            if not rel_id:
                break
            target = rel_map.get(rel_id)
            if not target:
                break
            target_path = posixpath.normpath(posixpath.join("xl", target))
            break
        if target_path is None:
            raise ValueError(f"{workbook_path} does not contain an `All_gene_pass_rates` sheet.")

        sheet_root = ET.fromstring(archive.read(target_path))
        sheet_rows = sheet_root.findall("main:sheetData/main:row", XLSX_MAIN_NS)
        if not sheet_rows:
            raise ValueError(f"{workbook_path} has an empty `All_gene_pass_rates` sheet.")

        def parse_row(row_element: ET.Element) -> List[str]:
            values: Dict[int, str] = {}
            max_index = -1
            for cell in row_element.findall("main:c", XLSX_MAIN_NS):
                ref = cell.attrib.get("r", "")
                column_index = excel_column_index(ref)
                max_index = max(max_index, column_index)
                cell_type = cell.attrib.get("t")
                if cell_type == "inlineStr":
                    inline_node = cell.find("main:is", XLSX_MAIN_NS)
                    values[column_index] = "".join(inline_node.itertext()) if inline_node is not None else ""
                    continue
                raw_value = cell.findtext("main:v", default="", namespaces=XLSX_MAIN_NS)
                if cell_type == "s":
                    try:
                        values[column_index] = shared_strings[int(raw_value)]
                    except (IndexError, ValueError):
                        values[column_index] = ""
                else:
                    values[column_index] = raw_value
            if max_index < 0:
                return []
            return [values.get(index, "") for index in range(max_index + 1)]

        headers = [normalize_gene_key(value) for value in parse_row(sheet_rows[0])]
        try:
            gene_column = headers.index("Gene")
        except ValueError as exc:
            raise ValueError(
                f"{workbook_path} `All_gene_pass_rates` is missing the `Gene` column."
            ) from exc

        if dataset_label:
            requested_column = f"{dataset_label}_Absolute_conversion_fraction"
            if requested_column not in headers:
                raise ValueError(
                    f"{workbook_path} does not contain the column `{requested_column}`."
                )
            ratio_columns = [headers.index(requested_column)]
        else:
            ratio_columns = [
                index
                for index, header in enumerate(headers)
                if header.endswith("_Absolute_conversion_fraction")
            ]
            if not ratio_columns:
                raise ValueError(
                    f"{workbook_path} `All_gene_pass_rates` has no `*_Absolute_conversion_fraction` columns."
                )

        allowed_genes: set[str] = set()
        for row_element in sheet_rows[1:]:
            row = parse_row(row_element)
            gene_name = normalize_gene_key(row[gene_column] if gene_column < len(row) else "")
            if not gene_name:
                continue
            for column_index in ratio_columns:
                if column_index >= len(row):
                    continue
                raw_value = row[column_index]
                if raw_value in ("", None):
                    continue
                try:
                    ratio_value = float(raw_value)
                except ValueError:
                    continue
                if ratio_value >= minimum_ratio:
                    allowed_genes.add(gene_name)
                    break
        return allowed_genes


def filter_transcript_models(
    transcript_models: Dict[str, TranscriptModel],
    allowed_gene_names: set[str] | None,
    explicit_gene_names: set[str] | None,
) -> Dict[str, TranscriptModel]:
    if not allowed_gene_names and not explicit_gene_names:
        return transcript_models

    filtered: Dict[str, TranscriptModel] = {}
    for transcript_id, model in transcript_models.items():
        gene_name = normalize_gene_key(model.gene_name)
        gene_id = normalize_gene_key(model.gene_id)
        transcript_key = normalize_gene_key(model.transcript_id)
        if allowed_gene_names:
            if gene_name not in allowed_gene_names and gene_id not in allowed_gene_names and transcript_key not in allowed_gene_names:
                continue
        if explicit_gene_names:
            if gene_name not in explicit_gene_names and gene_id not in explicit_gene_names and transcript_key not in explicit_gene_names:
                continue
        filtered[transcript_id] = model
    return filtered


def collect_coding_transcript_ids(gtf_path: Path, progress_interval: int) -> set[str]:
    coding_transcripts: set[str] = set()
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if transcript_id:
                coding_transcripts.add(transcript_id)
            if progress_interval > 0 and line_number % progress_interval == 0:
                print(
                    f"[gtf-pass1] lines={line_number:,} coding_transcripts={len(coding_transcripts):,}",
                    file=sys.stderr,
                    flush=True,
                )
    return coding_transcripts


def load_transcript_models(
    gtf_path: Path,
    progress_interval: int,
) -> Dict[str, TranscriptModel]:
    coding_tx_ids = collect_coding_transcript_ids(gtf_path, progress_interval)
    records: Dict[str, TranscriptRecord] = {}

    with gtf_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in {"transcript", "exon", "CDS"}:
                continue

            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id or transcript_id not in coding_tx_ids:
                continue

            record = records.setdefault(transcript_id, TranscriptRecord(transcript_id=transcript_id))
            chrom = fields[0]
            source = fields[1]
            start0 = int(fields[3]) - 1
            end0 = int(fields[4])
            strand = fields[6]
            record.chrom = record.chrom or chrom
            record.source = record.source or source
            record.strand = record.strand or strand
            record.gene_id = record.gene_id or normalize_label(attrs.get("gene_id"))
            record.gene_name = record.gene_name or normalize_label(attrs.get("gene"))

            if feature == "transcript":
                record.tx_start0 = start0 if record.tx_start0 is None else min(record.tx_start0, start0)
                record.tx_end0 = end0 if record.tx_end0 is None else max(record.tx_end0, end0)
            elif feature == "exon":
                record.exons.append((start0, end0))
            elif feature == "CDS":
                record.cds.append((start0, end0))

            if progress_interval > 0 and line_number % progress_interval == 0:
                print(
                    f"[gtf-pass2] lines={line_number:,} kept_transcripts={len(records):,}",
                    file=sys.stderr,
                    flush=True,
                )

    models: Dict[str, TranscriptModel] = {}
    for transcript_id, record in records.items():
        if not record.exons or not record.cds or not record.chrom:
            continue

        merged_exons = merge_intervals(record.exons)
        merged_cds = merge_intervals(record.cds)
        if not merged_exons or not merged_cds:
            continue

        tx_start0 = record.tx_start0 if record.tx_start0 is not None else merged_exons[0][0]
        tx_end0 = record.tx_end0 if record.tx_end0 is not None else merged_exons[-1][1]

        introns: List[Tuple[int, int]] = []
        for idx in range(len(merged_exons) - 1):
            intron_start0 = merged_exons[idx][1]
            intron_end0 = merged_exons[idx + 1][0]
            if intron_start0 < intron_end0:
                introns.append((intron_start0, intron_end0))

        cds_ordered = (
            sorted(merged_cds, reverse=True) if record.strand == "-" else sorted(merged_cds)
        )
        cds_chunks: List[CDSChunk] = []
        cds_offset0 = 0
        for start0, end0 in cds_ordered:
            cds_chunks.append(CDSChunk(start0=start0, end0=end0, cds_offset0=cds_offset0))
            cds_offset0 += end0 - start0

        cds_intervals = [(chunk.start0, chunk.end0) for chunk in cds_chunks]
        cds_intervals.sort()
        chunk_lookup = {(chunk.start0, chunk.end0): chunk for chunk in cds_chunks}
        ordered_chunks_by_genome = [chunk_lookup[interval] for interval in cds_intervals]

        cds_min0 = min(start0 for start0, _ in merged_cds)
        cds_max0 = max(end0 for _, end0 in merged_cds)

        models[transcript_id] = TranscriptModel(
            transcript_id=transcript_id,
            gene_id=record.gene_id or transcript_id,
            gene_name=record.gene_name or record.gene_id or transcript_id,
            chrom=record.chrom,
            strand=record.strand,
            source=record.source,
            tx_start0=tx_start0,
            tx_end0=tx_end0,
            exons=merged_exons,
            exon_starts=[start0 for start0, _ in merged_exons],
            introns=introns,
            intron_starts=[start0 for start0, _ in introns],
            cds_intervals=cds_intervals,
            cds_starts=[start0 for start0, _ in cds_intervals],
            cds_chunks=ordered_chunks_by_genome,
            cds_min0=cds_min0,
            cds_max0=cds_max0,
        )

    return models


def load_or_build_annotation_cache(
    gtf_path: Path,
    cache_path: Path | None,
    progress_interval: int,
) -> Tuple[Dict[str, TranscriptModel], Dict[str, Dict[int, List[str]]], bool]:
    if cache_path is not None and cache_path.exists():
        try:
            with cache_path.open("rb") as handle:
                payload = pickle.load(handle)
            stat = gtf_path.stat()
            if (
                payload.get("cache_version") == CACHE_VERSION
                and payload.get("gtf_path") == str(gtf_path.resolve())
                and payload.get("gtf_size") == stat.st_size
                and payload.get("gtf_mtime_ns") == stat.st_mtime_ns
            ):
                transcript_models = payload["transcript_models"]
                transcript_bucket_index = payload["transcript_bucket_index"]
                return transcript_models, transcript_bucket_index, True
        except Exception:
            pass

    transcript_models = load_transcript_models(gtf_path, progress_interval)
    transcript_bucket_index = build_transcript_bucket_index(transcript_models)

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        stat = gtf_path.stat()
        payload = {
            "cache_version": CACHE_VERSION,
            "gtf_path": str(gtf_path.resolve()),
            "gtf_size": stat.st_size,
            "gtf_mtime_ns": stat.st_mtime_ns,
            "transcript_models": transcript_models,
            "transcript_bucket_index": transcript_bucket_index,
        }
        with cache_path.open("wb") as handle:
            pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return transcript_models, transcript_bucket_index, False


def build_transcript_bucket_index(
    transcript_models: Dict[str, TranscriptModel],
) -> Dict[str, Dict[int, List[str]]]:
    index: DefaultDict[str, DefaultDict[int, List[str]]] = defaultdict(lambda: defaultdict(list))
    for transcript_id, model in transcript_models.items():
        first_bucket = model.tx_start0 // BUCKET_SIZE
        last_bucket = (model.tx_end0 - 1) // BUCKET_SIZE
        for bucket_id in range(first_bucket, last_bucket + 1):
            index[model.chrom][bucket_id].append(transcript_id)
    return {chrom: dict(bucket_map) for chrom, bucket_map in index.items()}


class PositionAnnotator:
    def __init__(
        self,
        transcript_models: Dict[str, TranscriptModel],
        transcript_bucket_index: Dict[str, Dict[int, List[str]]],
    ):
        self.transcript_models = transcript_models
        self.transcript_bucket_index = transcript_bucket_index
        self.bucket_caches: Dict[Tuple[str, int], Dict[int, int]] = {}
        self.c_position_buckets: Dict[Tuple[str, int], CPositionBucket] = {}
        self.last_bucket_key: Tuple[str, int] | None = None
        self.last_bucket_cache: Dict[int, int] | None = None
        self.last_chrom_buckets: Dict[int, List[str]] | None = None
        self.last_chrom: str | None = None
        self.last_c_bucket_key: Tuple[str, int] | None = None
        self.last_c_bucket: CPositionBucket | None = None

    @staticmethod
    def _encode_annotation(feature: str | None, codon_pos: int | None) -> int:
        if feature is None:
            return ANNOT_UNANNOTATED
        if feature == "5' UTR":
            return ANNOT_5UTR
        if feature == "3' UTR":
            return ANNOT_3UTR
        if feature == "intron":
            return ANNOT_INTRON
        if feature == "CDS":
            if codon_pos == 1:
                return ANNOT_CDS_1
            if codon_pos == 2:
                return ANNOT_CDS_2
            if codon_pos == 3:
                return ANNOT_CDS_3
            return ANNOT_CDS_AMBIGUOUS
        return ANNOT_AMBIGUOUS

    def block_has_any_annotation_candidates(self, chrom: str, start0: int, end0: int) -> bool:
        chrom_buckets = self._chrom_buckets(chrom)
        if not chrom_buckets:
            return False
        first_bucket = start0 // BUCKET_SIZE
        last_bucket = (end0 - 1) // BUCKET_SIZE
        for bucket_id in range(first_bucket, last_bucket + 1):
            if chrom_buckets.get(bucket_id):
                return True
        return False

    def _chrom_buckets(self, chrom: str) -> Dict[int, List[str]]:
        if chrom == self.last_chrom and self.last_chrom_buckets is not None:
            return self.last_chrom_buckets
        chrom_buckets = self.transcript_bucket_index.get(chrom, {})
        self.last_chrom = chrom
        self.last_chrom_buckets = chrom_buckets
        return chrom_buckets

    def annotate_code(self, chrom: str, pos0: int) -> int:
        bucket_id = pos0 // BUCKET_SIZE
        bucket_key = (chrom, bucket_id)
        if bucket_key == self.last_bucket_key and self.last_bucket_cache is not None:
            bucket_cache = self.last_bucket_cache
        else:
            bucket_cache = self.bucket_caches.get(bucket_key)
            if bucket_cache is None:
                bucket_cache = {}
                self.bucket_caches[bucket_key] = bucket_cache
                if len(self.bucket_caches) > 32:
                    self.bucket_caches.clear()
                    self.bucket_caches[bucket_key] = bucket_cache
            self.last_bucket_key = bucket_key
            self.last_bucket_cache = bucket_cache

        cached = bucket_cache.get(pos0)
        if cached is not None:
            return cached

        code = self._compute_annotation_code(chrom, pos0, bucket_id)
        bucket_cache[pos0] = code
        return code

    def _compute_annotation_code(self, chrom: str, pos0: int, bucket_id: int) -> int:
        candidate_ids = self._chrom_buckets(chrom).get(bucket_id, ())
        if not candidate_ids:
            return ANNOT_UNANNOTATED

        observed_feature: str | None = None
        observed_codon: int | None = None
        saw_any = False
        for transcript_id in candidate_ids:
            model = self.transcript_models[transcript_id]
            if pos0 < model.tx_start0 or pos0 >= model.tx_end0:
                continue
            feature, codon_pos = model.annotate_position(pos0)
            if feature is None:
                continue
            if not saw_any:
                observed_feature = feature
                observed_codon = codon_pos
                saw_any = True
                continue
            if feature != observed_feature:
                return ANNOT_AMBIGUOUS
            if feature == "CDS" and codon_pos != observed_codon:
                observed_codon = None

        if not saw_any:
            return ANNOT_UNANNOTATED
        return self._encode_annotation(observed_feature, observed_codon)

    def get_c_position_bucket(
        self,
        chrom: str,
        bucket_id: int,
        fasta: IndexedFasta,
    ) -> CPositionBucket:
        bucket_key = (chrom, bucket_id)
        if bucket_key == self.last_c_bucket_key and self.last_c_bucket is not None:
            return self.last_c_bucket

        cached = self.c_position_buckets.get(bucket_key)
        if cached is not None:
            self.last_c_bucket_key = bucket_key
            self.last_c_bucket = cached
            return cached

        bucket = self._build_c_position_bucket(chrom, bucket_id, fasta)
        self.c_position_buckets[bucket_key] = bucket
        if len(self.c_position_buckets) > 64:
            self.c_position_buckets.clear()
            self.c_position_buckets[bucket_key] = bucket
        self.last_c_bucket_key = bucket_key
        self.last_c_bucket = bucket
        return bucket

    def _build_c_position_bucket(
        self,
        chrom: str,
        bucket_id: int,
        fasta: IndexedFasta,
    ) -> CPositionBucket:
        candidate_ids = self._chrom_buckets(chrom).get(bucket_id, ())
        if not candidate_ids:
            return CPositionBucket(positions_by_code=[[] for _ in range(9)], code_by_position={})

        bucket_start0 = bucket_id * BUCKET_SIZE
        entry = fasta.index.get(chrom)
        if entry is None:
            raise KeyError(f"Reference sequence `{chrom}` not found in {fasta.fasta_path}")
        bucket_end0 = min(entry.length, bucket_start0 + BUCKET_SIZE)
        ref_seq = fasta.fetch(chrom, bucket_start0, bucket_end0)

        code_by_position: Dict[int, int] = {}
        for offset in iter_c_offsets(ref_seq):
            pos0 = bucket_start0 + offset
            code = self._compute_annotation_code(chrom, pos0, bucket_id)
            if code != ANNOT_UNANNOTATED:
                code_by_position[pos0] = code

        positions_by_code: List[List[int]] = [[] for _ in range(9)]
        for pos0, code in sorted(code_by_position.items()):
            positions_by_code[code].append(pos0)

        return CPositionBucket(
            positions_by_code=positions_by_code,
            code_by_position=code_by_position,
        )


def build_samtools_command(
    samtools: str,
    bam_path: Path,
    min_mapq: int,
    skip_flag: int,
    threads: int,
) -> List[str]:
    command = [samtools, "view", "-F", str(skip_flag)]
    if threads > 0:
        command.extend(["-@", str(threads)])
    if min_mapq > 0:
        command.extend(["-q", str(min_mapq)])
    command.append(str(bam_path))
    return command


def parse_nm_tag(fields: Sequence[str]) -> int | None:
    for tag_field in fields[11:]:
        if tag_field.startswith("nM:i:") or tag_field.startswith("NM:i:"):
            try:
                return int(tag_field[5:])
            except ValueError:
                return None
    return None


def iter_c_offsets(seq: str) -> Iterator[int]:
    offset = seq.find("C")
    while offset != -1:
        yield offset
        offset = seq.find("C", offset + 1)


def resolve_samtools_path(samtools: str) -> str:
    if os.path.sep in samtools:
        if Path(samtools).exists():
            return samtools
        raise FileNotFoundError(f"samtools not found at requested path: {samtools}")

    resolved = shutil.which(samtools)
    if resolved:
        return resolved

    common_locations = (
        "/opt/homebrew/bin/samtools",
        "/usr/local/bin/samtools",
        "/usr/bin/samtools",
    )
    for candidate in common_locations:
        if Path(candidate).exists():
            return candidate

    raise FileNotFoundError(
        "samtools was not found on PATH or in common install locations. "
        "Install it or pass --samtools /full/path/to/samtools."
    )


def estimate_primary_mapped_reads(samtools: str, bam_path: Path) -> int | None:
    command = [samtools, "idxstats", str(bam_path)]
    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            check=True,
        )
    except Exception:
        return None

    total = 0
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            total += int(fields[2])
    return total


def emit_progress(dataset: str, reads_seen: int, total_estimate: int | None, stats: DatasetStats) -> None:
    message = (
        f"[progress] dataset={dataset} reads={reads_seen:,} covered_C={stats.total_c_bases:,} "
        f"C_to_G={stats.total_c_to_g:,}"
    )
    if total_estimate:
        pct = min(100.0, (reads_seen / total_estimate) * 100.0)
        message += f" estimated_complete={pct:.2f}%"
    print(message, file=sys.stderr, flush=True)


def add_denominator_code_counts(stats: DatasetStats, code_counts: List[int], total_c_bases: int) -> None:
    cds_total = (
        code_counts[ANNOT_CDS_1]
        + code_counts[ANNOT_CDS_2]
        + code_counts[ANNOT_CDS_3]
        + code_counts[ANNOT_CDS_AMBIGUOUS]
    )
    annotated_or_ambiguous = (
        code_counts[ANNOT_5UTR]
        + cds_total
        + code_counts[ANNOT_3UTR]
        + code_counts[ANNOT_INTRON]
        + code_counts[ANNOT_AMBIGUOUS]
    )

    stats.total_c_bases += total_c_bases
    stats.feature_c_bases["5' UTR"] += code_counts[ANNOT_5UTR]
    stats.feature_c_bases["CDS"] += cds_total
    stats.feature_c_bases["3' UTR"] += code_counts[ANNOT_3UTR]
    stats.feature_c_bases["intron"] += code_counts[ANNOT_INTRON]
    stats.codon_c_bases[1] += code_counts[ANNOT_CDS_1]
    stats.codon_c_bases[2] += code_counts[ANNOT_CDS_2]
    stats.codon_c_bases[3] += code_counts[ANNOT_CDS_3]
    stats.ambiguous_codon_c_bases += code_counts[ANNOT_CDS_AMBIGUOUS]
    stats.ambiguous_feature_c_bases += code_counts[ANNOT_AMBIGUOUS]
    stats.unannotated_c_bases += max(0, total_c_bases - annotated_or_ambiguous)


def add_numerator_code(stats: DatasetStats, code: int) -> None:
    stats.total_c_to_g += 1
    if code == ANNOT_5UTR:
        stats.feature_c_to_g["5' UTR"] += 1
    elif code == ANNOT_3UTR:
        stats.feature_c_to_g["3' UTR"] += 1
    elif code == ANNOT_INTRON:
        stats.feature_c_to_g["intron"] += 1
    elif code == ANNOT_CDS_1:
        stats.feature_c_to_g["CDS"] += 1
        stats.codon_c_to_g[1] += 1
    elif code == ANNOT_CDS_2:
        stats.feature_c_to_g["CDS"] += 1
        stats.codon_c_to_g[2] += 1
    elif code == ANNOT_CDS_3:
        stats.feature_c_to_g["CDS"] += 1
        stats.codon_c_to_g[3] += 1
    elif code == ANNOT_CDS_AMBIGUOUS:
        stats.feature_c_to_g["CDS"] += 1
        stats.ambiguous_codon_c_to_g += 1
    elif code == ANNOT_AMBIGUOUS:
        stats.ambiguous_feature_c_to_g += 1
    else:
        stats.unannotated_c_to_g += 1


def analyze_bam(
    bam_path: Path,
    fasta: IndexedFasta,
    annotator: PositionAnnotator,
    samtools: str,
    samtools_threads: int,
    min_mapq: int,
    skip_flag: int,
    skip_zero_nm_reads: bool,
    max_reads: int | None,
    progress_interval: int,
) -> DatasetStats:
    stats = DatasetStats()
    total_estimate = estimate_primary_mapped_reads(samtools, bam_path)
    command = build_samtools_command(samtools, bam_path, min_mapq, skip_flag, samtools_threads)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )

    assert process.stdout is not None
    reached_limit = False
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            chrom = fields[2]
            pos1 = int(fields[3])
            cigar = fields[5]
            seq = fields[9].upper()
            if chrom == "*" or cigar == "*" or seq == "*":
                continue

            if max_reads is not None and stats.reads_seen >= max_reads:
                reached_limit = True
                break

            stats.reads_seen += 1
            nm = parse_nm_tag(fields)
            if skip_zero_nm_reads and nm == 0:
                stats.reads_skipped_nm_zero += 1
                continue

            try:
                nm_zero = nm == 0
                for read_start0, ref_start0, length in aligned_blocks(cigar, pos1 - 1):
                    block_end0 = ref_start0 + length
                    ref_seq = fasta.fetch(chrom, ref_start0, block_end0)
                    read_seq = seq[read_start0 : read_start0 + length]
                    read_has_only_bases = True
                    for nt in read_seq:
                        if nt not in BASES:
                            read_has_only_bases = False
                            break

                    if read_has_only_bases:
                        block_total_c = ref_seq.count("C")
                        code_counts = [0] * 9
                        first_bucket = ref_start0 // BUCKET_SIZE
                        last_bucket = (block_end0 - 1) // BUCKET_SIZE
                        for bucket_id in range(first_bucket, last_bucket + 1):
                            bucket = annotator.get_c_position_bucket(chrom, bucket_id, fasta)
                            sub_start0 = max(ref_start0, bucket_id * BUCKET_SIZE)
                            sub_end0 = min(block_end0, (bucket_id + 1) * BUCKET_SIZE)
                            for code in (
                                ANNOT_5UTR,
                                ANNOT_CDS_1,
                                ANNOT_CDS_2,
                                ANNOT_CDS_3,
                                ANNOT_CDS_AMBIGUOUS,
                                ANNOT_3UTR,
                                ANNOT_INTRON,
                                ANNOT_AMBIGUOUS,
                            ):
                                positions = bucket.positions_by_code[code]
                                if not positions:
                                    continue
                                lo = bisect_left(positions, sub_start0)
                                hi = bisect_left(positions, sub_end0)
                                if hi > lo:
                                    code_counts[code] += hi - lo
                        add_denominator_code_counts(stats, code_counts, block_total_c)
                    else:
                        block_total_c = 0
                        code_counts = [0] * 9
                        for offset in iter_c_offsets(ref_seq):
                            read_nt = read_seq[offset]
                            if read_nt not in BASES:
                                continue
                            block_total_c += 1
                            pos0 = ref_start0 + offset
                            bucket = annotator.get_c_position_bucket(chrom, pos0 // BUCKET_SIZE, fasta)
                            code = bucket.code_by_position.get(pos0, ANNOT_UNANNOTATED)
                            if code != ANNOT_UNANNOTATED:
                                code_counts[code] += 1
                        add_denominator_code_counts(stats, code_counts, block_total_c)

                    if nm_zero:
                        continue

                    for offset, read_nt in enumerate(read_seq):
                        if read_nt != "G" or ref_seq[offset] != "C":
                            continue
                        pos0 = ref_start0 + offset
                        bucket = annotator.get_c_position_bucket(chrom, pos0 // BUCKET_SIZE, fasta)
                        annotation_code = bucket.code_by_position.get(pos0, ANNOT_UNANNOTATED)
                        add_numerator_code(stats, annotation_code)
            except KeyError:
                stats.missing_reference_reads += 1
                continue
            except ValueError:
                continue

            if progress_interval > 0 and stats.reads_seen % progress_interval == 0:
                emit_progress(bam_path.name, stats.reads_seen, total_estimate, stats)
    finally:
        if process.stdout is not None:
            process.stdout.close()
        if reached_limit and process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if return_code != 0 and not reached_limit:
            raise RuntimeError(f"`{' '.join(command)}` failed with code {return_code}:\n{stderr}")

    return stats


def rate(numerator: int, denominator: int) -> float:
    return (numerator / denominator) if denominator else 0.0


def write_dataset_summary(
    output_path: Path,
    dataset_results: Dict[str, DatasetStats],
    combined: DatasetStats,
) -> None:
    headers = [
        "dataset",
        "reads_seen",
        "reads_skipped_nm_zero",
        "covered_reference_C_bases",
        "C_to_G_mismatches",
        "overall_C_to_G_rate",
        "annotated_feature_C_bases",
        "ambiguous_feature_C_bases",
        "ambiguous_feature_C_to_G",
        "unannotated_C_bases",
        "unannotated_C_to_G",
        "CDS_C_bases_with_resolved_codon",
        "CDS_C_to_G_with_resolved_codon",
        "CDS_C_bases_with_ambiguous_codon",
        "CDS_C_to_G_with_ambiguous_codon",
        "missing_reference_reads",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for dataset_name, stats in [*dataset_results.items(), ("ALL", combined)]:
            annotated_feature_c_bases = sum(stats.feature_c_bases.values())
            resolved_codon_c_bases = sum(stats.codon_c_bases.values())
            resolved_codon_c_to_g = sum(stats.codon_c_to_g.values())
            writer.writerow(
                    [
                        dataset_name,
                        stats.reads_seen,
                        stats.reads_skipped_nm_zero,
                        stats.total_c_bases,
                        stats.total_c_to_g,
                    f"{rate(stats.total_c_to_g, stats.total_c_bases):.8f}",
                    annotated_feature_c_bases,
                    stats.ambiguous_feature_c_bases,
                    stats.ambiguous_feature_c_to_g,
                    stats.unannotated_c_bases,
                    stats.unannotated_c_to_g,
                    resolved_codon_c_bases,
                    resolved_codon_c_to_g,
                    stats.ambiguous_codon_c_bases,
                    stats.ambiguous_codon_c_to_g,
                    stats.missing_reference_reads,
                ]
            )


def write_feature_summary(
    output_path: Path,
    dataset_results: Dict[str, DatasetStats],
    combined: DatasetStats,
) -> None:
    headers = [
        "dataset",
        "feature",
        "covered_reference_C_bases",
        "C_to_G_mismatches",
        "C_to_G_rate",
        "fraction_of_feature_annotated_C_to_G",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for dataset_name, stats in [*dataset_results.items(), ("ALL", combined)]:
            annotated_mismatches = sum(stats.feature_c_to_g.values())
            for feature in FEATURE_ORDER:
                covered_c = stats.feature_c_bases[feature]
                mismatches = stats.feature_c_to_g[feature]
                writer.writerow(
                    [
                        dataset_name,
                        feature,
                        covered_c,
                        mismatches,
                        f"{rate(mismatches, covered_c):.8f}",
                        f"{rate(mismatches, annotated_mismatches):.8f}",
                    ]
                )


def write_codon_summary(
    output_path: Path,
    dataset_results: Dict[str, DatasetStats],
    combined: DatasetStats,
) -> None:
    headers = [
        "dataset",
        "codon_position",
        "covered_reference_C_bases",
        "C_to_G_mismatches",
        "C_to_G_rate",
        "fraction_of_resolved_CDS_C_to_G",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for dataset_name, stats in [*dataset_results.items(), ("ALL", combined)]:
            resolved_cds_mismatches = sum(stats.codon_c_to_g.values())
            for codon_position in CODON_ORDER:
                covered_c = stats.codon_c_bases[codon_position]
                mismatches = stats.codon_c_to_g[codon_position]
                writer.writerow(
                    [
                        dataset_name,
                        codon_position,
                        covered_c,
                        mismatches,
                        f"{rate(mismatches, covered_c):.8f}",
                        f"{rate(mismatches, resolved_cds_mismatches):.8f}",
                    ]
                )


def write_metadata(
    output_path: Path,
    args: argparse.Namespace,
    bam_paths: Sequence[Path],
    transcript_models: Dict[str, TranscriptModel],
) -> None:
    per_chrom_counts: Counter = Counter(model.chrom for model in transcript_models.values())
    metadata = {
        "bam_paths": [str(path) for path in bam_paths],
        "gtf_path": str(args.gtf),
        "fasta_path": str(args.fasta),
        "samtools": args.samtools,
        "samtools_threads": args.samtools_threads,
        "min_mapq": args.min_mapq,
        "skip_flag": args.skip_flag,
        "skip_zero_nm_reads": args.skip_zero_nm_reads,
        "max_reads": args.max_reads,
        "progress_interval": args.progress_interval,
        "bucket_size": BUCKET_SIZE,
        "annotation_cache": str(args.annotation_cache) if args.annotation_cache else None,
        "notes": [
            "Both mates are counted as independent aligned reads.",
            "Only actually aligned/covered bases contribute to denominators.",
            "Counts are read-based, so the same genomic C can contribute multiple times if multiple reads cover it.",
            "Feature and codon summaries exclude ambiguous positions from denominators.",
            "Only coding transcripts with CDS entries in the GTF are used for feature/codon annotation.",
            "If skip_zero_nm_reads is enabled, reads with nM/NM equal to 0 are dropped entirely. This is faster but changes the denominator and should only be used for speed experiments.",
        ],
        "coding_transcript_count": len(transcript_models),
        "coding_transcripts_per_chromosome": dict(sorted(per_chrom_counts.items())),
    }
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)


def discover_bam_paths(bam_dir: Path, bam_glob: str) -> List[Path]:
    bam_paths = sorted(path for path in bam_dir.glob(bam_glob) if path.is_file())
    if not bam_paths:
        raise FileNotFoundError(f"No BAMs matched `{bam_glob}` in {bam_dir}")
    return bam_paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze C->G mismatches in RNA-seq BAMs by transcript feature "
            "(5' UTR, CDS, 3' UTR, intron) and CDS codon position."
        )
    )
    parser.add_argument(
        "--bam-dir",
        type=Path,
        default=Path("/Volumes/My Passport for Mac/ML012"),
        help="Directory containing coordinate-sorted BAM files [default: %(default)s]",
    )
    parser.add_argument(
        "--bam-glob",
        default="*_Aligned.sortedByCoord.out.bam",
        help="Glob used inside --bam-dir to find BAMs [default: %(default)s]",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        default=Path("/Volumes/My Passport for Mac/ML012/genomic.gtf"),
        help="Genome annotation GTF used to define CDS/UTR/intron [default: %(default)s]",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("/Volumes/My Passport for Mac/ML012/GCF_000001405.40_GRCh38.p14_genomic.fna"),
        help="Indexed genomic FASTA used to recover reference bases [default: %(default)s]",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("ML012_C2G_feature_codon"),
        help="Prefix for output files [default: %(default)s]",
    )
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="Path to samtools executable [default: %(default)s]",
    )
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=4,
        help="Threads to give samtools view via -@ [default: %(default)s]",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to pass to samtools view [default: %(default)s]",
    )
    parser.add_argument(
        "--skip-flag",
        type=lambda value: int(value, 0),
        default=DEFAULT_SKIP_FLAG,
        help=(
            "Bitmask of SAM flags to skip. Default removes unmapped, secondary, "
            "and supplementary alignments [default: %(default)s]"
        ),
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Optional per-BAM cap for quick tests",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=1_000_000,
        help="Emit stderr progress every N reads/lines while streaming [default: %(default)s]",
    )
    parser.add_argument(
        "--skip-zero-nm-reads",
        action="store_true",
        help=(
            "Experimental speed mode: drop reads whose nM/NM tag is 0 before any per-base work. "
            "This changes the denominator and is not appropriate for final rate estimates."
        ),
    )
    parser.add_argument(
        "--annotation-cache",
        type=Path,
        default=Path(".annotation_cache/genomic_gtf_coding_transcripts.pkl"),
        help="Pickle cache for parsed transcript annotation [default: %(default)s]",
    )
    parser.add_argument(
        "--gene-filter-workbook",
        type=Path,
        default=None,
        help="Workbook with `All_gene_pass_rates` used to keep only high-mismatch genes",
    )
    parser.add_argument(
        "--gene-filter-dataset",
        default=None,
        help="Specific dataset label in the workbook, for example `JM015M`; defaults to any dataset column",
    )
    parser.add_argument(
        "--gene-filter-min-ratio",
        type=float,
        default=0.05,
        help="Minimum workbook `Absolute_conversion_fraction` required to keep a gene [default: %(default)s]",
    )
    parser.add_argument(
        "--only-gene",
        action="append",
        default=[],
        help="Restrict analysis to this gene, gene_id, or transcript_id; may be passed multiple times",
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    if not (0.0 <= args.gene_filter_min_ratio <= 1.0):
        parser.error("--gene-filter-min-ratio must be between 0 and 1.")
    args.samtools = resolve_samtools_path(args.samtools)

    bam_paths = discover_bam_paths(args.bam_dir, args.bam_glob)
    print(f"[setup] discovered {len(bam_paths)} BAMs", file=sys.stderr, flush=True)
    print(f"[setup] using samtools at {args.samtools}", file=sys.stderr, flush=True)
    print(f"[setup] loading coding transcripts from {args.gtf}", file=sys.stderr, flush=True)
    transcript_models, transcript_bucket_index, used_cache = load_or_build_annotation_cache(
        args.gtf, args.annotation_cache, args.progress_interval
    )
    allowed_gene_names: set[str] | None = None
    if args.gene_filter_workbook is not None:
        allowed_gene_names = load_allowed_genes_from_workbook(
            args.gene_filter_workbook,
            args.gene_filter_dataset,
            args.gene_filter_min_ratio,
        )
        print(
            f"[setup] workbook filter kept {len(allowed_gene_names):,} genes at ratio >= {args.gene_filter_min_ratio:.4f}",
            file=sys.stderr,
            flush=True,
        )
    explicit_gene_names = {
        normalize_gene_key(gene_name) for gene_name in args.only_gene if normalize_gene_key(gene_name)
    }
    if explicit_gene_names:
        print(
            f"[setup] explicit gene restriction count={len(explicit_gene_names):,}",
            file=sys.stderr,
            flush=True,
        )
    transcript_models = filter_transcript_models(
        transcript_models,
        allowed_gene_names,
        explicit_gene_names if explicit_gene_names else None,
    )
    transcript_bucket_index = build_transcript_bucket_index(transcript_models)
    print(
        f"[setup] finalized {len(transcript_models):,} coding transcript models",
        file=sys.stderr,
        flush=True,
    )
    print(
        f"[setup] annotation cache {'hit' if used_cache else 'rebuilt'}",
        file=sys.stderr,
        flush=True,
    )
    annotator = PositionAnnotator(transcript_models, transcript_bucket_index)

    fasta = IndexedFasta(args.fasta)
    dataset_results: Dict[str, DatasetStats] = {}
    combined = DatasetStats()
    try:
        for bam_path in bam_paths:
            print(f"[bam] analyzing {bam_path.name}", file=sys.stderr, flush=True)
            stats = analyze_bam(
                bam_path=bam_path,
                fasta=fasta,
                annotator=annotator,
                samtools=args.samtools,
                samtools_threads=args.samtools_threads,
                min_mapq=args.min_mapq,
                skip_flag=args.skip_flag,
                skip_zero_nm_reads=args.skip_zero_nm_reads,
                max_reads=args.max_reads,
                progress_interval=args.progress_interval,
            )
            dataset_results[bam_path.stem] = stats
            combined.add(stats)
            print(
                f"[bam] finished {bam_path.name}: reads={stats.reads_seen:,} covered_C={stats.total_c_bases:,} C_to_G={stats.total_c_to_g:,}",
                file=sys.stderr,
                flush=True,
            )
    finally:
        fasta.close()

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    dataset_summary_path = Path(f"{output_prefix}.dataset_summary.tsv")
    feature_summary_path = Path(f"{output_prefix}.feature_summary.tsv")
    codon_summary_path = Path(f"{output_prefix}.codon_summary.tsv")
    metadata_path = Path(f"{output_prefix}.run_metadata.json")

    write_dataset_summary(dataset_summary_path, dataset_results, combined)
    write_feature_summary(feature_summary_path, dataset_results, combined)
    write_codon_summary(codon_summary_path, dataset_results, combined)
    write_metadata(metadata_path, args, bam_paths, transcript_models)

    print(f"dataset summary: {dataset_summary_path}", file=sys.stderr, flush=True)
    print(f"feature summary: {feature_summary_path}", file=sys.stderr, flush=True)
    print(f"codon summary: {codon_summary_path}", file=sys.stderr, flush=True)
    print(f"run metadata: {metadata_path}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
