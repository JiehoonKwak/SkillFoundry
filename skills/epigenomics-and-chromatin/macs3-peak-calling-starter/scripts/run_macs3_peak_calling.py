#!/usr/bin/env python3
"""Run a deterministic toy MACS3 peak-calling example."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MACS3 = ROOT / "slurm" / "envs" / "genomics" / "bin" / "macs3"


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_peak_xls(path: Path) -> tuple[int, int]:
    total_tags = 0
    fragment_length = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# total tags in treatment:"):
            total_tags = int(line.rsplit(":", 1)[1].strip())
        if line.startswith("# d ="):
            fragment_length = int(float(line.rsplit("=", 1)[1].strip()))
    return total_tags, fragment_length


def parse_narrowpeak(path: Path) -> list[dict[str, object]]:
    peaks: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            peaks.append(
                {
                    "chrom": row[0],
                    "start": int(row[1]),
                    "end": int(row[2]),
                    "name": row[3],
                    "score": int(row[4]),
                    "strand": row[5],
                    "signal_value": round(float(row[6]), 4),
                    "log10_pvalue": round(float(row[7]), 4),
                    "log10_qvalue": round(float(row[8]), 4),
                    "summit_offset": int(row[9]),
                }
            )
    return peaks


def summarize_callpeak(treatment_path: Path, workdir: Path) -> dict[str, object]:
    if not MACS3.exists():
        raise SystemExit("Missing MACS3 executable in slurm/envs/genomics.")
    workdir.mkdir(parents=True, exist_ok=True)
    outdir = workdir / "macs3"
    outdir.mkdir(parents=True, exist_ok=True)
    name = "toy"
    subprocess.run(
        [
            str(MACS3),
            "callpeak",
            "-t",
            str(treatment_path),
            "-f",
            "BED",
            "-g",
            "1000",
            "-q",
            "0.5",
            "--nomodel",
            "--extsize",
            "75",
            "--keep-dup",
            "all",
            "--name",
            name,
            "--outdir",
            str(outdir),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=240,
    )

    xls_path = outdir / f"{name}_peaks.xls"
    narrowpeak_path = outdir / f"{name}_peaks.narrowPeak"
    summits_path = outdir / f"{name}_summits.bed"
    peaks = parse_narrowpeak(narrowpeak_path)
    total_tags, fragment_length = parse_peak_xls(xls_path)
    payload = {
        "treatment_path": display_path(treatment_path),
        "peak_count": len(peaks),
        "treatment_tags": total_tags,
        "fragment_length": fragment_length,
        "outputs": {
            "narrowpeak": display_path(narrowpeak_path),
            "summits_bed": display_path(summits_path),
            "peaks_xls": display_path(xls_path),
        },
        "top_peak": peaks[0] if peaks else None,
    }
    return payload


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--treatment", type=Path, required=True, help="BED file of toy treatment intervals.")
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("scratch/epigenomics/macs3_peak_calling"),
        help="Work directory for MACS3 outputs.",
    )
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    treatment_path = args.treatment if args.treatment.is_absolute() else ROOT / args.treatment
    workdir = args.workdir if args.workdir.is_absolute() else ROOT / args.workdir
    out_path = args.summary_out if args.summary_out is None or args.summary_out.is_absolute() else ROOT / args.summary_out
    payload = summarize_callpeak(treatment_path.resolve(), workdir.resolve())
    write_json(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
