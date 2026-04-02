#!/usr/bin/env python3
"""Render a conservative Slurm batch script."""

from __future__ import annotations

import argparse
from pathlib import Path


def render_script(job_name: str, command: str, partition: str | None, time_limit: str, mem: str, cpus: int) -> str:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --output=slurm/logs/%x-%j.out",
        "#SBATCH --error=slurm/logs/%x-%j.err",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={mem}",
    ]
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "echo \"hostname: $(hostname)\"",
            "echo \"started: $(date --iso-8601=seconds)\"",
            command,
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--command", required=True, help="Command to run inside the batch job")
    parser.add_argument("--job-name", default="smoke", help="Slurm job name")
    parser.add_argument("--partition", default=None, help="Optional partition name")
    parser.add_argument("--time", dest="time_limit", default="00:05:00", help="SBATCH time limit")
    parser.add_argument("--mem", default="1G", help="SBATCH memory request")
    parser.add_argument("--cpus", type=int, default=1, help="SBATCH cpus-per-task")
    parser.add_argument("--out", type=Path, default=None, help="Optional output path")
    args = parser.parse_args()

    script = render_script(args.job_name, args.command, args.partition, args.time_limit, args.mem, args.cpus)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(script, encoding="utf-8")
    else:
        print(script, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
