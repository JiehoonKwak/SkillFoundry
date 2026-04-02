#!/usr/bin/env python3
"""Summarize a SAM/BAM/CRAM file with pysam and optionally emit an indexed BAM copy."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pysam


def infer_input_mode(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".sam":
        return "r"
    if suffix == ".bam":
        return "rb"
    if suffix == ".cram":
        return "rc"
    raise SystemExit(f"Unsupported input extension for {path.name}; expected .sam, .bam, or .cram")


def summarize_alignment_file(input_path: Path, bam_out: Path | None) -> dict[str, object]:
    mode = infer_input_mode(input_path)
    reference_mapped_counts: Counter[str] = Counter()
    record_count = 0
    mapped_count = 0
    unmapped_count = 0
    paired_count = 0
    proper_pair_count = 0
    secondary_count = 0
    supplementary_count = 0
    query_length_sum = 0

    with pysam.AlignmentFile(str(input_path), mode) as reader:
        reference_names = list(reader.references)
        writer = None
        if bam_out is not None:
            bam_out.parent.mkdir(parents=True, exist_ok=True)
            writer = pysam.AlignmentFile(str(bam_out), "wb", template=reader)
        try:
            for record in reader.fetch(until_eof=True):
                record_count += 1
                query_length_sum += record.query_length or 0
                if record.is_paired:
                    paired_count += 1
                if record.is_proper_pair:
                    proper_pair_count += 1
                if record.is_secondary:
                    secondary_count += 1
                if record.is_supplementary:
                    supplementary_count += 1
                if record.is_unmapped:
                    unmapped_count += 1
                else:
                    mapped_count += 1
                    reference_name = reader.get_reference_name(record.reference_id)
                    reference_mapped_counts[reference_name] += 1
                if writer is not None:
                    writer.write(record)
        finally:
            if writer is not None:
                writer.close()

    bai_path = None
    if bam_out is not None:
        pysam.index(str(bam_out))
        bai_path = Path(f"{bam_out}.bai")

    mean_query_length = round(query_length_sum / record_count, 6) if record_count else 0.0
    return {
        "input_path": str(input_path),
        "input_format": mode,
        "reference_count": len(reference_names),
        "reference_names": reference_names,
        "record_count": record_count,
        "mapped_count": mapped_count,
        "unmapped_count": unmapped_count,
        "paired_count": paired_count,
        "proper_pair_count": proper_pair_count,
        "secondary_count": secondary_count,
        "supplementary_count": supplementary_count,
        "mean_query_length": mean_query_length,
        "reference_mapped_counts": dict(sorted(reference_mapped_counts.items())),
        "bam_written": bam_out is not None,
        "bam_path": str(bam_out) if bam_out is not None else None,
        "bam_index_path": str(bai_path) if bai_path is not None else None,
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Input SAM, BAM, or CRAM path.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--bam-out", type=Path, default=None, help="Optional BAM output path to write and index.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    payload = summarize_alignment_file(args.input, args.bam_out)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
