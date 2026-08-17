#!/usr/bin/env python3
"""
How to run:

python3 calculate_gene_mismatch_halflife.py \
  --sam "/path/to/*F_trim2.SAM" \
  --gtf "/path/to/refseq.gtf" \
  --half-life-csv "/path/to/input_values.csv" \
  --offset 12 \
  --min-cds-reads 20 \
  --out "/path/to/output.xlsx"

Required flags:
  --sam            One or more SAM files or glob patterns.
  --gtf            GTF file with transcript/CDS annotations.
  --half-life-csv  One or more CSV files with headers "Gene Name" and "Half Life".

Optional flags:
  --offset         Offset from the read 5' end. Default: 12
  --min-cds-reads  Minimum CDS-eligible reads required to report a gene. Default: 0
  --out            Output path. Use .xlsx for multi-sheet workbook output or .csv for one SAM input.
"""

import argparse
import csv
import glob
import math
import os
import re
import sys
from bisect import bisect_right
from collections import defaultdict

from openpyxl import Workbook
from scipy.stats import pearsonr, spearmanr


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_gtf_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, val = part.split(" ", 1)
            attrs[key] = val.strip().strip('"')
        elif "=" in part:
            key, val = part.split("=", 1)
            attrs[key] = val.strip().strip('"')
    return attrs


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def build_search_arrays(intervals):
    starts = [start for start, _ in intervals]
    ends = [end for _, end in intervals]
    return starts, ends


def in_intervals(pos, starts, ends):
    if not starts:
        return False
    idx = bisect_right(starts, pos) - 1
    return idx >= 0 and pos <= ends[idx]


def ref_aligned_length(cigar):
    if cigar == "*" or not cigar:
        return 0
    length = 0
    for count, op in _CIGAR_RE.findall(cigar):
        if op in ("M", "D", "N", "=", "X"):
            length += int(count)
    return length


def five_prime_transcript_pos(flag, pos, cigar):
    if flag & 0x10:
        return pos + ref_aligned_length(cigar) - 1
    return pos


def apply_offset_from_5prime(flag, five_prime, offset):
    if flag & 0x10:
        return five_prime - offset
    return five_prime + offset


def normalize_id(value):
    if not value:
        return None
    return value.split(".", 1)[0]


def add_alias(alias_map, alias, value):
    if alias:
        alias_map[alias].add(value)
        normalized = normalize_id(alias)
        if normalized and normalized != alias:
            alias_map[normalized].add(value)


def parse_optional_int_tag(fields, tag):
    prefix = f"{tag}:"
    for field in fields[11:]:
        if field.startswith(prefix):
            parts = field.split(":", 2)
            if len(parts) == 3:
                try:
                    return int(parts[2])
                except ValueError:
                    return None
    return None


def has_mismatch(fields):
    value = parse_optional_int_tag(fields, "NM")
    return value is not None and value > 0


def build_transcript_cds_index(gtf_path, gene_attr="gene_id", transcript_attr="transcript_id"):
    exons_by_tx = defaultdict(list)
    cds_by_tx = defaultdict(list)
    strand_by_tx = {}
    gene_name_by_gene_id = {}
    gene_name_by_tx = {}
    tx_to_gene_id = {}
    alias_to_txs = defaultdict(set)
    alias_to_gene_names = defaultdict(set)

    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqname, _, feature, start, end, _, strand, _, attrs_str = parts
            if feature not in {"gene", "transcript", "exon", "CDS"}:
                continue
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            attrs = parse_gtf_attributes(attrs_str)
            gene_id = attrs.get(gene_attr) or attrs.get("gene_id") or attrs.get("gene")
            tx_id = attrs.get(transcript_attr) or attrs.get("transcript_id") or attrs.get("transcript")
            gene_name = attrs.get("gene_name") or seqname or gene_id
            tx_name = attrs.get("transcript_name")

            if gene_id and gene_name:
                gene_name_by_gene_id.setdefault(gene_id, gene_name)
                add_alias(alias_to_gene_names, gene_id, gene_name)
                add_alias(alias_to_gene_names, gene_name, gene_name)
                add_alias(alias_to_gene_names, seqname, gene_name)

            if feature == "gene":
                continue

            if not tx_id:
                continue

            strand_by_tx[tx_id] = strand
            if gene_id:
                tx_to_gene_id[tx_id] = gene_id
            if gene_name:
                gene_name_by_tx[tx_id] = gene_name

            add_alias(alias_to_txs, tx_id, tx_id)
            add_alias(alias_to_txs, seqname, tx_id)
            if tx_name:
                add_alias(alias_to_txs, tx_name, tx_id)
            if gene_id:
                add_alias(alias_to_txs, gene_id, tx_id)
            if gene_name:
                add_alias(alias_to_txs, gene_name, tx_id)

            if feature == "exon":
                exons_by_tx[tx_id].append((start, end))
            elif feature == "CDS":
                cds_by_tx[tx_id].append((start, end))

    cds_search_by_tx = {}
    gene_name_by_tx_final = {}

    for tx_id in set(exons_by_tx) | set(cds_by_tx):
        exons = exons_by_tx.get(tx_id)
        cds_parts = cds_by_tx.get(tx_id)
        if not exons or not cds_parts:
            continue
        strand = strand_by_tx.get(tx_id, "+")
        exons = sorted(exons, reverse=(strand == "-"))

        exon_tx_coords = []
        tx_cursor = 1
        for exon_start, exon_end in exons:
            exon_len = exon_end - exon_start + 1
            tx_start = tx_cursor
            tx_end = tx_cursor + exon_len - 1
            exon_tx_coords.append((exon_start, exon_end, tx_start, tx_end))
            tx_cursor += exon_len

        cds_tx_intervals = []
        for cds_start, cds_end in cds_parts:
            for exon_start, exon_end, tx_start, _ in exon_tx_coords:
                ov_start = max(cds_start, exon_start)
                ov_end = min(cds_end, exon_end)
                if ov_start > ov_end:
                    continue
                if strand == "+":
                    mapped_start = tx_start + (ov_start - exon_start)
                    mapped_end = tx_start + (ov_end - exon_start)
                else:
                    mapped_start = tx_start + (exon_end - ov_end)
                    mapped_end = tx_start + (exon_end - ov_start)
                if mapped_start > mapped_end:
                    mapped_start, mapped_end = mapped_end, mapped_start
                cds_tx_intervals.append((mapped_start, mapped_end))

        merged = merge_intervals(cds_tx_intervals)
        if not merged:
            continue
        cds_search_by_tx[tx_id] = build_search_arrays(merged)
        gene_id = tx_to_gene_id.get(tx_id)
        gene_name = gene_name_by_tx.get(tx_id) or gene_name_by_gene_id.get(gene_id) or gene_id or tx_id
        gene_name_by_tx_final[tx_id] = gene_name

    return cds_search_by_tx, gene_name_by_tx_final, alias_to_txs, alias_to_gene_names


def resolve_gene_name(rname, gene_name_by_tx, alias_to_txs, alias_to_gene_names):
    tx_candidates = alias_to_txs.get(rname)
    if not tx_candidates:
        tx_candidates = alias_to_txs.get(normalize_id(rname), set())
    if len(tx_candidates) == 1:
        tx_id = next(iter(tx_candidates))
        return tx_id, gene_name_by_tx.get(tx_id)

    gene_candidates = alias_to_gene_names.get(rname)
    if not gene_candidates:
        gene_candidates = alias_to_gene_names.get(normalize_id(rname), set())
    if len(gene_candidates) == 1:
        return None, next(iter(gene_candidates))

    if len(tx_candidates) > 1:
        for tx_id in sorted(tx_candidates):
            gene_name = gene_name_by_tx.get(tx_id)
            if gene_name:
                return tx_id, gene_name
    return None, None


def load_average_half_lives(paths):
    values_by_gene = defaultdict(list)
    for path in paths:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                gene_name = (row.get("Gene Name") or "").strip()
                half_life = (row.get("Half Life") or "").strip()
                if not gene_name or not half_life:
                    continue
                try:
                    values_by_gene[gene_name].append(float(half_life))
                except ValueError:
                    continue
    return {
        gene_name: sum(values) / len(values)
        for gene_name, values in values_by_gene.items()
        if values
    }


def compute_correlations(base_rows):
    ratios = []
    log_half_lives = []
    for _, ratio, half_life in base_rows:
        if half_life == "" or half_life is None:
            continue
        if half_life <= 0:
            continue
        ratios.append(ratio)
        log_half_lives.append(math.log10(half_life))

    if len(ratios) < 2:
        return "", ""

    pearson_value = pearsonr(ratios, log_half_lives).statistic
    spearman_value = spearmanr(ratios, log_half_lives).statistic
    return pearson_value, spearman_value


def build_output_rows(total_by_gene, mismatch_by_gene, avg_half_life, min_cds_reads):
    base_rows = []
    for gene_name in sorted(total_by_gene):
        total_reads = total_by_gene[gene_name]
        if total_reads < min_cds_reads:
            continue
        mismatch_value = mismatch_by_gene.get(gene_name, 0)
        ratio = mismatch_value / total_reads if total_reads else 0.0
        base_rows.append((gene_name, ratio, avg_half_life.get(gene_name, "")))

    pearson_value, spearman_value = compute_correlations(base_rows)
    rows = []
    for gene_name, ratio, half_life in base_rows:
        log_half_life = math.log10(half_life) if half_life not in ("", None) and half_life > 0 else ""
        rows.append((gene_name, ratio, half_life, log_half_life, spearman_value, pearson_value))
    return rows


def collect_gene_stats(sam_path, cds_search_by_tx, gene_name_by_tx, alias_to_txs, alias_to_gene_names,
                       offset):
    total_by_gene = defaultdict(int)
    mismatch_by_gene = defaultdict(int)
    unresolved = set()
    processed = 0
    kept = 0

    with open(sam_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            processed += 1
            try:
                flag = int(fields[1])
                pos = int(fields[3])
            except ValueError:
                continue

            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue

            rname = fields[2]
            cigar = fields[5]
            if rname == "*" or cigar == "*":
                continue

            tx_id, gene_name = resolve_gene_name(rname, gene_name_by_tx, alias_to_txs, alias_to_gene_names)
            if tx_id is None or gene_name is None:
                unresolved.add(rname)
                continue

            starts, ends = cds_search_by_tx.get(tx_id, (None, None))
            if starts is None:
                unresolved.add(rname)
                continue

            five_prime = five_prime_transcript_pos(flag, pos, cigar)
            site = apply_offset_from_5prime(flag, five_prime, offset)
            if not in_intervals(site, starts, ends):
                continue

            kept += 1
            total_by_gene[gene_name] += 1

            if has_mismatch(fields):
                mismatch_by_gene[gene_name] += 1

    return total_by_gene, mismatch_by_gene, processed, kept, unresolved


def write_csv_output(out_path, rows):
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "gene_name",
            "mismatch_ratio",
            "average_half_life",
            "log10_average_half_life",
            "spearman_correlation",
            "pearson_correlation",
        ])
        for idx, (gene_name, ratio, half_life, log_half_life, spearman_value, pearson_value) in enumerate(rows):
            writer.writerow([
                gene_name,
                f"{ratio:.10f}",
                half_life,
                log_half_life,
                spearman_value if idx == 0 else "",
                pearson_value if idx == 0 else "",
            ])


def safe_sheet_name(name, used_names):
    cleaned = re.sub(r"[\[\]\*\:/\\\?]", "_", name).strip("'")
    cleaned = cleaned or "Sheet"
    cleaned = cleaned[:31]
    candidate = cleaned
    suffix = 1
    while candidate in used_names:
        suffix_text = f"_{suffix}"
        candidate = f"{cleaned[:31 - len(suffix_text)]}{suffix_text}"
        suffix += 1
    used_names.add(candidate)
    return candidate


def write_excel_output(out_path, sheets):
    wb = Workbook()
    default_ws = wb.active
    wb.remove(default_ws)
    used_names = set()

    for sheet_base_name, rows in sheets:
        ws = wb.create_sheet(title=safe_sheet_name(sheet_base_name, used_names))
        ws.append([
            "gene_name",
            "mismatch_ratio",
            "average_half_life",
            "log10_average_half_life",
            "spearman_correlation",
            "pearson_correlation",
        ])
        for idx, (gene_name, ratio, half_life, log_half_life, spearman_value, pearson_value) in enumerate(rows):
            ws.append([
                gene_name,
                ratio,
                half_life,
                log_half_life,
                spearman_value if idx == 0 else "",
                pearson_value if idx == 0 else "",
            ])

    wb.save(out_path)


def expand_sam_paths(patterns):
    paths = []
    seen = set()
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            for path in matches:
                if path not in seen:
                    paths.append(path)
                    seen.add(path)
        elif os.path.exists(pattern):
            if pattern not in seen:
                paths.append(pattern)
                seen.add(pattern)
    return paths


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate per-gene mismatch ratios from a transcriptome-aligned SAM, restricted to reads "
            "whose 5' end plus offset lands inside the CDS, and join averaged half-life values."
        )
    )
    parser.add_argument(
        "--sam",
        required=True,
        nargs="+",
        help="One or more SAM files or glob patterns, for example '*F_trim2.SAM'",
    )
    parser.add_argument("--gtf", required=True, help="GTF with transcript/CDS annotations")
    parser.add_argument(
        "--half-life-csv",
        dest="half_life_csvs",
        required=True,
        nargs="+",
        help="One or more half-life replicate CSV files.",
    )
    parser.add_argument("--offset", type=int, default=12, help="Offset from read 5' end in nt")
    parser.add_argument(
        "--min-cds-reads",
        type=int,
        default=0,
        help="Optional cutoff: only report genes with at least this many CDS-eligible reads",
    )
    parser.add_argument(
        "--out",
        default="gene_mismatch_halflife.xlsx",
        help="Output path. Use .xlsx for multi-sheet workbook, or .csv for a single SAM file.",
    )
    args = parser.parse_args()

    sam_paths = expand_sam_paths(args.sam)
    if not sam_paths:
        print("No SAM files matched the provided --sam inputs.", file=sys.stderr)
        return 1

    half_life_csvs = args.half_life_csvs

    missing_inputs = [path for path in [args.gtf, *half_life_csvs] if not os.path.exists(path)]
    if missing_inputs:
        for path in missing_inputs:
            print(f"Missing input: {path}", file=sys.stderr)
        return 1

    cds_search_by_tx, gene_name_by_tx, alias_to_txs, alias_to_gene_names = build_transcript_cds_index(args.gtf)
    avg_half_life = load_average_half_lives(half_life_csvs)

    sheets = []
    for sam_path in sam_paths:
        total_by_gene, mismatch_by_gene, processed, kept, unresolved = collect_gene_stats(
            sam_path=sam_path,
            cds_search_by_tx=cds_search_by_tx,
            gene_name_by_tx=gene_name_by_tx,
            alias_to_txs=alias_to_txs,
            alias_to_gene_names=alias_to_gene_names,
            offset=args.offset,
        )
        rows = build_output_rows(total_by_gene, mismatch_by_gene, avg_half_life, args.min_cds_reads)
        sheets.append((os.path.splitext(os.path.basename(sam_path))[0], rows))

        print(f"{os.path.basename(sam_path)}", file=sys.stderr)
        print(f"  Processed alignments: {processed}", file=sys.stderr)
        print(f"  CDS-eligible alignments counted: {kept}", file=sys.stderr)
        print(f"  Genes written: {len(rows)}", file=sys.stderr)
        if unresolved:
            preview = ", ".join(sorted(unresolved)[:10])
            print(
                f"  Unresolved reference names skipped: {len(unresolved)}"
                + (f" ({preview})" if preview else ""),
                file=sys.stderr,
            )

    if args.out.lower().endswith(".csv"):
        if len(sheets) != 1:
            print("CSV output supports exactly one SAM input. Use .xlsx for multiple SAM files.", file=sys.stderr)
            return 1
        write_csv_output(args.out, sheets[0][1])
    else:
        write_excel_output(args.out, sheets)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
