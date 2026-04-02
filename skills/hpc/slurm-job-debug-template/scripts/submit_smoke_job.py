#!/usr/bin/env python3
"""Submit a tiny Slurm smoke job and capture its accounting summary."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from render_sbatch import render_script


ROOT = Path(__file__).resolve().parents[4]


def default_command(sleep_seconds: int) -> str:
    return (
        "echo hostname=$(hostname)\n"
        "echo started=$(date --iso-8601=seconds)\n"
        f"sleep {sleep_seconds}\n"
        "echo finished=$(date --iso-8601=seconds)"
    )


def parse_sacct_table(text: str) -> dict[str, str]:
    lines = [line for line in text.strip().splitlines() if line.strip()]
    if len(lines) < 2:
        return {}
    headers = lines[0].split("|")
    values = lines[1].split("|")
    return dict(zip(headers, values, strict=True))


def fetch_accounting(job_id: str, attempts: int = 5) -> dict[str, str]:
    for attempt in range(attempts):
        completed = subprocess.run(
            [
                "sacct",
                "-X",
                "-j",
                job_id,
                "--format=JobIDRaw,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList",
                "-P",
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        parsed = parse_sacct_table(completed.stdout)
        if parsed:
            return parsed
        if attempt + 1 < attempts:
            time.sleep(1)
    raise RuntimeError(f"sacct did not return accounting for job {job_id}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition", default="cpu", help="Slurm partition for the smoke job")
    parser.add_argument("--job-name", default="slurm-smoke", help="Slurm job name")
    parser.add_argument("--time", dest="time_limit", default="00:02:00", help="SBATCH time limit")
    parser.add_argument("--mem", default="512M", help="SBATCH memory request")
    parser.add_argument("--cpus", type=int, default=1, help="SBATCH cpus-per-task")
    parser.add_argument("--sleep", type=int, default=2, help="Sleep time inside the smoke job")
    parser.add_argument("--command", default=None, help="Optional custom command body")
    parser.add_argument("--script-out", type=Path, default=None, help="Optional path to save the rendered sbatch script")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON report path")
    args = parser.parse_args()

    (ROOT / "slurm" / "jobs").mkdir(parents=True, exist_ok=True)
    (ROOT / "slurm" / "logs").mkdir(parents=True, exist_ok=True)
    (ROOT / "slurm" / "reports").mkdir(parents=True, exist_ok=True)

    script_path = args.script_out or ROOT / "slurm" / "jobs" / f"{args.job_name}.sbatch"
    command = args.command or default_command(args.sleep)
    script_text = render_script(
        job_name=args.job_name,
        command=command,
        partition=args.partition,
        time_limit=args.time_limit,
        mem=args.mem,
        cpus=args.cpus,
    )
    script_path.write_text(script_text, encoding="utf-8")

    submit = subprocess.run(
        ["sbatch", "--parsable", "--wait", str(script_path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
    )
    job_id = submit.stdout.strip().split(";", 1)[0]
    accounting = fetch_accounting(job_id)
    stdout_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}.out"
    stderr_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}.err"

    payload = {
        "job_id": job_id,
        "job_name": args.job_name,
        "partition": args.partition,
        "script_path": str(script_path.relative_to(ROOT)),
        "stdout_path": str(stdout_path.relative_to(ROOT)),
        "stderr_path": str(stderr_path.relative_to(ROOT)),
        "accounting": accounting,
        "stdout": stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else "",
        "stderr": stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else "",
    }

    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
