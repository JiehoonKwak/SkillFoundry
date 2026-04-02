#!/usr/bin/env python3
"""Filter a small VCF with bcftools and summarize the retained variants."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
BCFTOOLS = ROOT / "slurm" / "envs" / "bcftools" / "bin" / "bcftools"
BCFTOOLS_LIB = ROOT / "slurm" / "envs" / "bcftools" / "lib"


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def tool_env() -> dict[str, str]:
    env = os.environ.copy()
    ld_library_path = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = f"{BCFTOOLS_LIB}:{ld_library_path}" if ld_library_path else str(BCFTOOLS_LIB)
    return env


def record_ids_from_vcf(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    ids: list[str] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            ids.append(line.split("\t", 3)[2])
    return ids


def summarize_filter(input_path: Path, filtered_vcf_out: Path, include_expr: str) -> dict[str, object]:
    if not BCFTOOLS.exists():
        raise SystemExit("Missing bcftools executable in slurm/envs/bcftools.")

    filtered_vcf_out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(BCFTOOLS),
            "view",
            "--include",
            include_expr,
            "--output-type",
            "z",
            "--output-file",
            str(filtered_vcf_out),
            str(input_path),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=240,
        env=tool_env(),
    )
    subprocess.run(
        [str(BCFTOOLS), "index", "--tbi", str(filtered_vcf_out)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=240,
        env=tool_env(),
    )

    input_ids = record_ids_from_vcf(input_path)
    kept_ids = record_ids_from_vcf(filtered_vcf_out)
    index_path = Path(f"{filtered_vcf_out}.tbi")
    return {
        "input_path": display_path(input_path),
        "filter_expression": include_expr,
        "input_record_count": len(input_ids),
        "passing_record_count": len(kept_ids),
        "kept_ids": kept_ids,
        "filtered_vcf_path": display_path(filtered_vcf_out),
        "filtered_vcf_index_path": display_path(index_path),
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
    parser.add_argument("--input", type=Path, required=True, help="Input VCF path.")
    parser.add_argument(
        "--include-expr",
        default="QUAL>=50 && INFO/DP>=15",
        help="bcftools include expression used with `bcftools view --include`.",
    )
    parser.add_argument(
        "--filtered-vcf-out",
        type=Path,
        default=Path("scratch/genomics/toy_variants.filtered.vcf.gz"),
        help="Output compressed filtered VCF path.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    input_path = args.input if args.input.is_absolute() else ROOT / args.input
    filtered_vcf_out = args.filtered_vcf_out if args.filtered_vcf_out.is_absolute() else ROOT / args.filtered_vcf_out
    out_path = args.out if args.out is None or args.out.is_absolute() else ROOT / args.out

    if not input_path.exists():
        raise SystemExit(f"Input VCF not found: {input_path}")

    payload = summarize_filter(input_path.resolve(), filtered_vcf_out.resolve(), args.include_expr)
    write_json(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
