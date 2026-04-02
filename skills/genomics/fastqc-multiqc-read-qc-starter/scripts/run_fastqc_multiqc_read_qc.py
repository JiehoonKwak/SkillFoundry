#!/usr/bin/env python3
"""Run FastQC and MultiQC on a small FASTQ input and summarize the result."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import zipfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
ENV_BIN = ROOT / "slurm" / "envs" / "genomics" / "bin"
FASTQC = ENV_BIN / "fastqc"
MULTIQC = ENV_BIN / "multiqc"


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def build_runtime_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    env.setdefault("JAVA_HOME", str(ENV_BIN.parent))
    return env


def require_executable(path: Path, name: str) -> None:
    if not path.exists():
        raise SystemExit(
            f"Missing {name} at {path}. Install the genomics prefix before running this skill."
        )


def run_command(args: list[str], env: dict[str, str]) -> None:
    subprocess.run(
        args,
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
        timeout=240,
    )


def extract_fastqc_members(zip_path: Path) -> tuple[list[str], dict[str, str], list[tuple[str, str, str]]]:
    with zipfile.ZipFile(zip_path) as handle:
        data_member = next(name for name in handle.namelist() if name.endswith("/fastqc_data.txt"))
        summary_member = next(name for name in handle.namelist() if name.endswith("/summary.txt"))
        fastqc_data = handle.read(data_member).decode("utf-8")
        summary_rows = [
            tuple(line.split("\t"))
            for line in handle.read(summary_member).decode("utf-8").splitlines()
            if line.strip()
        ]
    basic_stats: dict[str, str] = {}
    in_basic_stats = False
    for line in fastqc_data.splitlines():
        if line.startswith(">>Basic Statistics"):
            in_basic_stats = True
            continue
        if in_basic_stats and line.startswith(">>END_MODULE"):
            break
        if in_basic_stats and line and not line.startswith("#"):
            key, value = line.split("\t", 1)
            basic_stats[key] = value
    return fastqc_data.splitlines(), basic_stats, summary_rows


def parse_multiqc_general_stats(path: Path) -> tuple[int, list[str], dict[str, str] | None]:
    if not path.exists():
        return 0, [], None
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    if not rows:
        return 0, [], None
    sample_names = [row["Sample"] for row in rows if row.get("Sample")]
    first_row = dict(rows[0])
    return len(rows), sample_names, first_row


def summarize_read_qc(input_path: Path, workdir: Path) -> dict:
    require_executable(FASTQC, "FastQC")
    require_executable(MULTIQC, "MultiQC")
    env = build_runtime_env()
    if shutil.which("java", path=env["PATH"]) is None:
        raise SystemExit("FastQC runtime is not usable because java is missing from the genomics prefix PATH.")

    sample_name = input_path.stem
    fastqc_out = workdir / "fastqc"
    multiqc_out = workdir / "multiqc"
    fastqc_out.mkdir(parents=True, exist_ok=True)
    multiqc_out.mkdir(parents=True, exist_ok=True)

    run_command([str(FASTQC), "-o", str(fastqc_out), str(input_path)], env=env)
    run_command([str(MULTIQC), "-o", str(multiqc_out), str(fastqc_out)], env=env)

    fastqc_zip = fastqc_out / f"{sample_name}_fastqc.zip"
    fastqc_html = fastqc_out / f"{sample_name}_fastqc.html"
    if not fastqc_zip.exists():
        raise FileNotFoundError(f"Expected FastQC zip output at {fastqc_zip}")
    _, basic_stats, summary_rows = extract_fastqc_members(fastqc_zip)
    module_counts = dict(sorted(Counter(status for status, _, _ in summary_rows).items()))

    general_stats_path = multiqc_out / "multiqc_data" / "multiqc_general_stats.txt"
    multiqc_sample_count, multiqc_samples, multiqc_first_row = parse_multiqc_general_stats(general_stats_path)
    payload = {
        "sample_name": sample_name,
        "input_fastq": display_path(input_path),
        "total_sequences": int(basic_stats["Total Sequences"]),
        "total_bases": basic_stats["Total Bases"],
        "sequence_length": basic_stats["Sequence length"],
        "gc_percent": int(basic_stats["%GC"]),
        "fastqc_module_status_counts": module_counts,
        "multiqc_sample_count": multiqc_sample_count,
        "multiqc_samples": multiqc_samples,
        "reports": {
            "fastqc_zip": display_path(fastqc_zip),
            "fastqc_html": display_path(fastqc_html),
            "multiqc_html": display_path(multiqc_out / "multiqc_report.html"),
        },
        "multiqc_general_stats_row": multiqc_first_row,
    }
    return payload


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="FASTQ file to analyze.")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("scratch/genomics/fastqc_multiqc_work"),
        help="Working directory for FastQC and MultiQC outputs.",
    )
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    input_path = args.input
    if not input_path.is_absolute():
        input_path = ROOT / input_path
    payload = summarize_read_qc(input_path.resolve(), args.workdir if args.workdir.is_absolute() else ROOT / args.workdir)
    write_json(payload, args.summary_out if args.summary_out is None or args.summary_out.is_absolute() else ROOT / args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
