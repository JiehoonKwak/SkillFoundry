#!/usr/bin/env python3
"""Submit a tiny Slurm job, monitor it with squeue, and summarize final accounting."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]


def render_probe_script(job_name: str, partition: str, sleep_seconds: int, cpus: int, mem: str, time_limit: str) -> str:
    return "\n".join(
        [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --partition={partition}",
            f"#SBATCH --time={time_limit}",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --mem={mem}",
            f"#SBATCH --output={ROOT / 'slurm' / 'logs' / (job_name + '-%j.out')}",
            f"#SBATCH --error={ROOT / 'slurm' / 'logs' / (job_name + '-%j.err')}",
            "echo hostname=$(hostname)",
            "echo started=$(date --iso-8601=seconds)",
            f"sleep {sleep_seconds}",
            "echo finished=$(date --iso-8601=seconds)",
        ]
    ) + "\n"


def parse_pipe_row(text: str, fields: list[str]) -> dict[str, str] | None:
    row = text.strip().split("|")
    if len(row) != len(fields):
        return None
    return dict(zip(fields, row, strict=True))


def submit_probe(job_name: str, partition: str, sleep_seconds: int, cpus: int, mem: str, time_limit: str) -> tuple[str, Path]:
    script_path = ROOT / "slurm" / "jobs" / f"{job_name}-monitoring.sbatch"
    script_path.write_text(
        render_probe_script(
            job_name=job_name,
            partition=partition,
            sleep_seconds=sleep_seconds,
            cpus=cpus,
            mem=mem,
            time_limit=time_limit,
        ),
        encoding="utf-8",
    )
    completed = subprocess.run(
        ["sbatch", "--parsable", str(script_path)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    job_id = completed.stdout.strip().split(";", 1)[0]
    return job_id, script_path


def collect_squeue_snapshots(job_id: str, poll_interval: float, max_polls: int) -> list[dict[str, str]]:
    snapshots: list[dict[str, str]] = []
    fields = ["JobID", "State", "Partition", "TimeUsed", "Nodes", "NodeOrReason"]
    for _ in range(max_polls):
        completed = subprocess.run(
            ["squeue", "-h", "-j", job_id, "-o", "%i|%T|%P|%M|%D|%R"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        text = completed.stdout.strip()
        if text:
            parsed = parse_pipe_row(text, fields)
            if parsed is not None:
                snapshots.append(parsed)
        else:
            if snapshots:
                break
        time.sleep(poll_interval)
    return snapshots


ACTIVE_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "RESIZING", "SUSPENDED"}


def fetch_sacct(job_id: str, attempts: int = 30) -> dict[str, str]:
    fields = ["JobIDRaw", "JobName", "Partition", "State", "ExitCode", "Elapsed", "AllocCPUS", "NodeList"]
    for attempt in range(attempts):
        completed = subprocess.run(
            [
                "sacct",
                "-X",
                "-n",
                "-P",
                "-j",
                job_id,
                "--format=JobIDRaw,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,NodeList",
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        if lines:
            parsed = parse_pipe_row(lines[0], fields)
            if parsed is not None and parsed["State"] not in ACTIVE_STATES:
                return parsed
        if attempt + 1 < attempts:
            time.sleep(1)
    raise RuntimeError(f"sacct did not return terminal accounting for job {job_id}")


def try_seff(job_id: str) -> str | None:
    if shutil.which("seff") is None:
        return None
    completed = subprocess.run(
        ["seff", job_id],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    lines = []
    for line in (completed.stdout + completed.stderr).splitlines():
        line = line.strip()
        if not line or "Use of uninitialized value $user" in line:
            continue
        lines.append(line)
    text = "\n".join(lines).strip()
    return text or None


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition", default="cpu", help="Slurm partition for the probe job.")
    parser.add_argument("--job-name", default="slurm-monitoring", help="Slurm job name.")
    parser.add_argument("--sleep", type=int, default=2, help="Sleep time inside the probe job.")
    parser.add_argument("--cpus", type=int, default=1, help="SBATCH cpus-per-task.")
    parser.add_argument("--mem", default="256M", help="SBATCH memory request.")
    parser.add_argument("--time", dest="time_limit", default="00:02:00", help="SBATCH time limit.")
    parser.add_argument("--poll-interval", type=float, default=1.0, help="Polling interval for squeue.")
    parser.add_argument("--max-polls", type=int, default=10, help="Maximum number of squeue polls.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON report path.")
    args = parser.parse_args()

    if args.max_polls <= 0:
        raise SystemExit("max-polls must be a positive integer.")
    if args.poll_interval <= 0:
        raise SystemExit("poll-interval must be positive.")

    (ROOT / "slurm" / "jobs").mkdir(parents=True, exist_ok=True)
    (ROOT / "slurm" / "logs").mkdir(parents=True, exist_ok=True)

    job_id, script_path = submit_probe(
        job_name=args.job_name,
        partition=args.partition,
        sleep_seconds=args.sleep,
        cpus=args.cpus,
        mem=args.mem,
        time_limit=args.time_limit,
    )
    queue_snapshots = collect_squeue_snapshots(job_id=job_id, poll_interval=args.poll_interval, max_polls=args.max_polls)
    accounting = fetch_sacct(job_id)
    stdout_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}.out"
    stderr_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}.err"
    payload = {
        "job_id": job_id,
        "job_name": args.job_name,
        "partition": args.partition,
        "script_path": str(script_path.relative_to(ROOT)),
        "queue_snapshot_count": len(queue_snapshots),
        "queue_states_seen": sorted({snapshot["State"] for snapshot in queue_snapshots}),
        "queue_snapshots": queue_snapshots,
        "accounting": accounting,
        "seff_summary": try_seff(job_id),
        "stdout_path": str(stdout_path.relative_to(ROOT)),
        "stderr_path": str(stderr_path.relative_to(ROOT)),
        "stdout": stdout_path.read_text(encoding="utf-8") if stdout_path.exists() else "",
        "stderr": stderr_path.read_text(encoding="utf-8") if stderr_path.exists() else "",
    }
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
