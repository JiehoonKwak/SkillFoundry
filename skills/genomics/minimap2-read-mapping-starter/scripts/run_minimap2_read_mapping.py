#!/usr/bin/env python3
"""Run a deterministic tiny minimap2 short-read mapping workflow."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "examples" / "toy_reference.fa"
DEFAULT_READS = ROOT / "examples" / "toy_reads.fastq"


def tool_path(name: str) -> Path:
    path = Path(__file__).resolve().parents[4] / "slurm" / "envs" / "genomics" / "bin" / name
    if not path.exists():
        raise SystemExit(f"Required tool not found: {path}")
    return path


def count_fastq_reads(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        line_count = sum(1 for _ in handle)
    if line_count % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {path}")
    return line_count // 4


def build_summary(reference: Path, reads: Path, bam_out: Path) -> dict[str, object]:
    minimap2 = tool_path("minimap2")
    samtools = tool_path("samtools")
    bam_out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="minimap2-", dir=bam_out.parent) as tmp_dir:
        tmp_root = Path(tmp_dir)
        sam_path = tmp_root / "mapped.sam"
        with sam_path.open("w", encoding="utf-8") as handle:
            subprocess.run(
                [str(minimap2), "-a", "-x", "sr", str(reference), str(reads)],
                check=True,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
            )
        subprocess.run(
            [str(samtools), "sort", "-o", str(bam_out), str(sam_path)],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        subprocess.run(
            [str(samtools), "index", str(bam_out)],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )

    alignments: list[dict[str, object]] = []
    mapped_count = 0
    total_mapq = 0
    reference_names: set[str] = set()
    with pysam.AlignmentFile(str(bam_out), "rb") as bam_file:
        reference_lengths = {
            name: int(length)
            for name, length in zip(bam_file.references, bam_file.lengths)
        }
        for record in bam_file.fetch(until_eof=True):
            if not record.is_unmapped:
                mapped_count += 1
                total_mapq += record.mapping_quality
                if record.reference_name:
                    reference_names.add(record.reference_name)
            alignments.append(
                {
                    "read_name": record.query_name,
                    "mapped": not record.is_unmapped,
                    "reference_name": record.reference_name,
                    "start": None if record.is_unmapped else int(record.reference_start),
                    "mapq": int(record.mapping_quality),
                    "cigar": record.cigarstring,
                }
            )

    read_count = count_fastq_reads(reads)
    return {
        "reference_path": str(reference),
        "reads_path": str(reads),
        "bam_path": str(bam_out),
        "bai_path": str(bam_out) + ".bai",
        "read_count": read_count,
        "mapped_count": mapped_count,
        "unmapped_count": read_count - mapped_count,
        "mean_mapq": round(total_mapq / mapped_count, 6) if mapped_count else 0.0,
        "reference_names": sorted(reference_names),
        "reference_lengths": reference_lengths,
        "alignments": alignments,
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
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE, help="Reference FASTA path.")
    parser.add_argument("--reads", type=Path, default=DEFAULT_READS, help="FASTQ reads path.")
    parser.add_argument("--bam-out", type=Path, required=True, help="Sorted BAM output path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    if not args.reference.exists():
        raise SystemExit(f"Reference FASTA not found: {args.reference}")
    if not args.reads.exists():
        raise SystemExit(f"FASTQ reads not found: {args.reads}")
    if args.bam_out.exists():
        args.bam_out.unlink()
    bai_path = Path(str(args.bam_out) + ".bai")
    if bai_path.exists():
        bai_path.unlink()

    payload = build_summary(args.reference, args.reads, args.bam_out)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
