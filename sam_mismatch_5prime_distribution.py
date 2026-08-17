#!/usr/bin/env python3
"""
Build 5' mismatch-distribution workbooks from SAM files.

Example:
    python3 sam_mismatch_5prime_distribution.py \
        -d /path/to/sam_directory \
        -o sam_mismatch_5prime_distribution.xlsx \
        --cds-gtf refseq.gtf \
        --offset 12

High-level flow:
1. Load CDS bounds from GTF.
2. Keep only reads whose shifted 5' position (with user offset) falls in CDS.
3. Parse mismatches and count position-wise distributions from the read 5' end.
4. Write:
   - 4 global sheets (from_A/from_C/from_G/from_T),
   - one per-dataset detailed sheet with from->to composition by position.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import zipfile
from collections import Counter, defaultdict
from typing import DefaultDict, Dict, List, Tuple
from xml.sax.saxutils import escape


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MD_RE = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")
BASES = ("A", "C", "G", "T")
REF_SPAN_OPS = ("M", "D", "N", "=", "X")


def parse_tags(fields: List[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for tag_field in fields[11:]:
        parts = tag_field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_gtf_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key_value = part.split(" ", 1)
        if len(key_value) != 2:
            continue
        key, value = key_value
        attrs[key] = value.strip().strip('"')
    return attrs


def normalize_ref_name(name: str) -> str:
    base = name.strip().split("|", 1)[0]
    if base.endswith("_transcript"):
        return base[: -len("_transcript")]
    return base


def load_cds_ranges(gtf_path: str) -> Dict[str, Tuple[int, int]]:
    """Build transcript/gene -> [min CDS start, max CDS end] from GTF CDS rows."""
    cds_ranges: Dict[str, List[int]] = {}
    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")

            names = set()
            if gene_id:
                names.add(normalize_ref_name(gene_id))
            if transcript_id:
                names.add(normalize_ref_name(transcript_id))

            for name in names:
                if not name:
                    continue
                if name not in cds_ranges:
                    cds_ranges[name] = [start, end]
                else:
                    cds_ranges[name][0] = min(cds_ranges[name][0], start)
                    cds_ranges[name][1] = max(cds_ranges[name][1], end)

    return {k: (v[0], v[1]) for k, v in cds_ranges.items()}


def reference_span(cigar: str) -> int:
    span = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        if op in REF_SPAN_OPS:
            span += int(length_str)
    return span


def shifted_five_prime_position(pos: int, cigar: str, is_reverse: bool, offset: int) -> int:
    """
    Compute the offset-adjusted 5' genomic position for CDS gating.

    Forward read:  shifted = 5prime + offset
    Reverse read:  shifted = 5prime - offset
    """
    span = reference_span(cigar)
    if span <= 0:
        return -1
    five_prime = (pos + span - 1) if is_reverse else pos
    return (five_prime - offset) if is_reverse else (five_prime + offset)


def build_ref_to_read_map(cigar: str) -> Dict[int, int]:
    ref_to_read: Dict[int, int] = {}
    read_pos = 0
    ref_pos = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str)
        if op in ("M", "=", "X"):
            for i in range(length):
                ref_to_read[ref_pos + i] = read_pos + i
            read_pos += length
            ref_pos += length
        elif op in ("I", "S"):
            read_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue
    return ref_to_read


def extract_mismatch_events(seq: str, cigar: str, md: str) -> List[Tuple[str, str, int]]:
    """Return (ref_base, read_base, read_index_0_based) mismatch events."""
    ref_to_read = build_ref_to_read_map(cigar)
    mismatches: List[Tuple[str, str, int]] = []
    ref_offset = 0
    seq = seq.upper()

    for token in MD_RE.findall(md.upper()):
        if token[0].isdigit():
            ref_offset += int(token)
        elif token.startswith("^"):
            ref_offset += len(token) - 1
        else:
            ref_base = token
            read_idx = ref_to_read.get(ref_offset)
            if read_idx is not None and 0 <= read_idx < len(seq):
                read_base = seq[read_idx]
                if ref_base in BASES and read_base in BASES and ref_base != read_base:
                    mismatches.append((ref_base, read_base, read_idx))
            ref_offset += 1

    return mismatches


def analyze_dataset(
    path: str, cds_ranges: Dict[str, Tuple[int, int]], offset: int
) -> Dict[str, Dict]:
    """
    For each from-base, count mismatch events at 5' positions.
    Returns:
      - position_counts: from_base -> Counter(position_1_based_from_5prime -> count)
      - transition_counts: (from_base, to_base) -> Counter(position_1_based_from_5prime -> count)
    """
    position_counts: DefaultDict[str, Counter] = defaultdict(Counter)
    transition_counts: DefaultDict[Tuple[str, str], Counter] = defaultdict(Counter)

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            rname = normalize_ref_name(fields[2])
            pos = int(fields[3])
            cigar = fields[5]
            seq = fields[9]
            if (flag & 4) != 0 or cigar == "*" or seq == "*":
                continue

            # CDS inclusion filter: only analyze reads whose shifted 5' position is in CDS.
            cds = cds_ranges.get(rname)
            is_reverse = (flag & 16) != 0
            if cds is None:
                continue
            shifted_pos = shifted_five_prime_position(pos, cigar, is_reverse, offset)
            if shifted_pos < cds[0] or shifted_pos > cds[1]:
                continue

            tags = parse_tags(fields)
            md = tags.get("MD")
            if not md:
                continue

            read_len = len(seq)
            # p1 is closest to the read 5' end (strand-aware).
            for ref_base, read_base, read_idx in extract_mismatch_events(seq, cigar, md):
                pos5 = (read_len - read_idx) if is_reverse else (read_idx + 1)
                position_counts[ref_base][pos5] += 1
                transition_counts[(ref_base, read_base)][pos5] += 1

    return {
        "position_counts": position_counts,
        "transition_counts": transition_counts,
    }


def excel_col_name(index_1_based: int) -> str:
    result = ""
    index = index_1_based
    while index > 0:
        index, rem = divmod(index - 1, 26)
        result = chr(ord("A") + rem) + result
    return result


def xml_cell(value: object) -> str:
    if isinstance(value, (int, float)):
        return f'<c t="n"><v>{value}</v></c>'
    return f'<c t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'


def worksheet_xml(headers: List[str], rows: List[List[object]]) -> str:
    last_col = excel_col_name(len(headers))
    last_row = len(rows) + 1
    sheet_rows = []
    sheet_rows.append("<row r=\"1\">" + "".join(xml_cell(h) for h in headers) + "</row>")
    for i, row in enumerate(rows, start=2):
        sheet_rows.append(f"<row r=\"{i}\">" + "".join(xml_cell(v) for v in row) + "</row>")
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="A1:{last_col}{last_row}"/>'
        f"<sheetData>{''.join(sheet_rows)}</sheetData>"
        "</worksheet>"
    )


def write_xlsx_multi(output_path: str, sheets: List[Tuple[str, List[str], List[List[object]]]]) -> None:
    """Write a minimal multi-sheet XLSX (no external Excel deps)."""
    workbook_sheet_xml = []
    workbook_rels_xml = []
    content_type_overrides = [
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
    ]

    for i, (name, _headers, _rows) in enumerate(sheets, start=1):
        workbook_sheet_xml.append(
            f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
        )
        workbook_rels_xml.append(
            f'<Relationship Id="rId{i}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{i}.xml"/>'
        )
        content_type_overrides.append(
            f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{''.join(workbook_sheet_xml)}</sheets>"
        "</workbook>"
    )

    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f"{''.join(content_type_overrides)}"
        "</Types>"
    )

    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )

    workbook_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{''.join(workbook_rels_xml)}"
        "</Relationships>"
    )

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", root_rels_xml)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        for i, (_name, headers, rows) in enumerate(sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{i}.xml", worksheet_xml(headers, rows))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a workbook of 5' mismatch distributions (percent) across SAM datasets."
    )
    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        help="Directory containing SAM files (default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="sam_mismatch_5prime_distribution_all.xlsx",
        help="Output Excel filename (default: sam_mismatch_5prime_distribution_all.xlsx)",
    )
    parser.add_argument(
        "--cds-gtf",
        default=None,
        help="GTF file containing CDS features (prompted if omitted)",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=None,
        help="Offset from read 5' position used for CDS inclusion test (prompted if omitted)",
    )
    args = parser.parse_args()

    sam_paths = sorted(
        glob.glob(os.path.join(args.directory, "*.SAM"))
        + glob.glob(os.path.join(args.directory, "*.sam"))
    )
    if not sam_paths:
        raise SystemExit(f"No SAM files found in: {os.path.abspath(args.directory)}")
    cds_gtf = args.cds_gtf
    if cds_gtf is None:
        cds_gtf = input("Enter path to CDS GTF file: ").strip()
    if not cds_gtf:
        raise SystemExit("No CDS GTF path provided.")
    if not os.path.exists(cds_gtf):
        raise SystemExit(f"CDS GTF file not found: {os.path.abspath(cds_gtf)}")

    offset = args.offset
    if offset is None:
        offset = int(input("Enter offset for CDS inclusion filtering: ").strip())

    cds_ranges = load_cds_ranges(cds_gtf)
    if not cds_ranges:
        raise SystemExit(f"No CDS records found in: {os.path.abspath(cds_gtf)}")

    dataset_names = [os.path.basename(p) for p in sam_paths]
    dataset_data = {
        os.path.basename(p): analyze_dataset(p, cds_ranges, offset) for p in sam_paths
    }
    max_pos = max(
        (
            max(
                dataset_data[ds]["position_counts"].get(from_base, Counter()).keys(),
                default=0,
            )
            for ds in dataset_names
            for from_base in BASES
        ),
        default=0,
    )

    # Global sheets: position distribution for each mutated-from base.
    sheets: List[Tuple[str, List[str], List[List[object]]]] = []
    for from_base in BASES:
        headers = ["position_5prime"] + [f"{ds}_percent" for ds in dataset_names]

        rows: List[List[object]] = []
        for pos in range(1, max_pos + 1):
            row: List[object] = [pos]
            for ds in dataset_names:
                pos_counter = dataset_data[ds]["position_counts"].get(from_base, Counter())
                total_from_base = sum(pos_counter.values())
                if total_from_base == 0:
                    percent = 0.0
                else:
                    percent = (pos_counter.get(pos, 0) / total_from_base) * 100.0
                row.append(percent)
            rows.append(row)

        sheets.append((f"from_{from_base}", headers, rows))

    # Detailed sheets: within each dataset and position, split each from-base by to-base.
    transition_pairs = [
        (from_base, to_base)
        for from_base in BASES
        for to_base in BASES
        if from_base != to_base
    ]
    for ds in dataset_names:
        ds_stem = os.path.splitext(ds)[0]
        sheet_name = f"detail_{ds_stem}"[:31]
        headers = ["position_5prime"] + [f"{f}_to_{t}_percent" for f, t in transition_pairs]
        rows: List[List[object]] = []
        for pos in range(1, max_pos + 1):
            row: List[object] = [pos]
            for from_base, to_base in transition_pairs:
                from_total_at_pos = dataset_data[ds]["position_counts"].get(
                    from_base, Counter()
                ).get(pos, 0)
                trans_at_pos = dataset_data[ds]["transition_counts"].get(
                    (from_base, to_base), Counter()
                ).get(pos, 0)
                if from_total_at_pos == 0:
                    pct = 0.0
                else:
                    pct = (trans_at_pos / from_total_at_pos) * 100.0
                row.append(pct)
            rows.append(row)
        sheets.append((sheet_name, headers, rows))

    write_xlsx_multi(args.output, sheets)
    print(f"Wrote 5' distribution workbook to: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
