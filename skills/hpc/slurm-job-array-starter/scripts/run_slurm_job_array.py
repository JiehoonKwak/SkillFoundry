#!/usr/bin/env python3
"""Submit a tiny Slurm job array and summarize per-task accounting."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
ACTIVE_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "RESIZING", "SUSPENDED"}


def expand_array_spec(spec: str) -> list[int]:
    tasks: set[int] = set()
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" not in chunk:
            tasks.add(int(chunk))
            continue
        range_part, _, step_part = chunk.partition(":")
        start_text, end_text = range_part.split("-", 1)
        start = int(start_text)
        end = int(end_text)
        step = int(step_part) if step_part else 1
        if step <= 0:
            raise ValueError("Array step must be a positive integer.")
        if end < start:
            raise ValueError("Array range end must be >= start.")
        tasks.update(range(start, end + 1, step))
    if not tasks:
        raise ValueError("Array spec must contain at least one task.")
    return sorted(tasks)


def render_array_script(job_name: str, partition: str, array_spec: str, sleep_seconds: int, cpus: int, mem: str, time_limit: str) -> str:
    return "\n".join(
        [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --partition={partition}",
            f"#SBATCH --time={time_limit}",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --mem={mem}",
            f"#SBATCH --array={array_spec}",
            f"#SBATCH --output={ROOT / 'slurm' / 'logs' / (job_name + '-%A_%a.out')}",
            f"#SBATCH --error={ROOT / 'slurm' / 'logs' / (job_name + '-%A_%a.err')}",
            "echo array_job_id=$SLURM_ARRAY_JOB_ID",
            "echo array_task_id=$SLURM_ARRAY_TASK_ID",
            "echo started=$(date --iso-8601=seconds)",
            f"sleep {sleep_seconds}",
            "echo finished=$(date --iso-8601=seconds)",
        ]
    ) + "\n"


def parse_pipe_row(text: str, fields: list[str]) -> dict[str, str] | None:
    parts = text.strip().split("|")
    if len(parts) != len(fields):
        return None
    return dict(zip(fields, parts, strict=True))


def submit_array_job(job_name: str, partition: str, array_spec: str, sleep_seconds: int, cpus: int, mem: str, time_limit: str) -> tuple[str, Path]:
    script_path = ROOT / "slurm" / "jobs" / f"{job_name}-array.sbatch"
    script_path.write_text(
        render_array_script(
            job_name=job_name,
            partition=partition,
            array_spec=array_spec,
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


def collect_squeue_snapshots(job_id: str, poll_interval: float, max_polls: int) -> list[list[dict[str, str]]]:
    fields = ["JobID", "State", "TimeUsed", "NodeOrReason"]
    snapshots: list[list[dict[str, str]]] = []
    for _ in range(max_polls):
        completed = subprocess.run(
            ["squeue", "-h", "-j", job_id, "-o", "%i|%T|%M|%R"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        rows = []
        for line in completed.stdout.splitlines():
            parsed = parse_pipe_row(line, fields)
            if parsed is not None:
                rows.append(parsed)
        if rows:
            snapshots.append(rows)
        elif snapshots:
            break
        time.sleep(poll_interval)
    return snapshots


def fetch_task_accounting(job_id: str, task_id: int, attempts: int = 120) -> dict[str, str]:
    fields = ["JobIDRaw", "JobName", "State", "ExitCode", "Elapsed"]
    query_id = f"{job_id}_{task_id}"
    for attempt in range(attempts):
        completed = subprocess.run(
            ["sacct", "-n", "-P", "-j", query_id, "--format=JobIDRaw,JobName,State,ExitCode,Elapsed"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        rows = []
        for line in completed.stdout.splitlines():
            if not line.strip() or ".batch" in line or ".extern" in line:
                continue
            parsed = parse_pipe_row(line, fields)
            if parsed is not None:
                rows.append(parsed)
        if rows and rows[0]["State"] not in ACTIVE_STATES:
            return rows[0]
        if attempt + 1 < attempts:
            time.sleep(1)
    raise RuntimeError(f"sacct did not return terminal accounting for array task {query_id}")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition", default="cpu", help="Slurm partition for the array job.")
    parser.add_argument("--job-name", default="slurm-array", help="Array job name.")
    parser.add_argument("--array-spec", default="0-1", help="Slurm array specification, for example 0-1 or 0-4:2.")
    parser.add_argument("--sleep", type=int, default=1, help="Sleep time inside each array task.")
    parser.add_argument("--cpus", type=int, default=1, help="SBATCH cpus-per-task.")
    parser.add_argument("--mem", default="256M", help="SBATCH memory request.")
    parser.add_argument("--time", dest="time_limit", default="00:02:00", help="SBATCH time limit.")
    parser.add_argument("--poll-interval", type=float, default=1.0, help="Queue polling interval.")
    parser.add_argument("--max-polls", type=int, default=20, help="Maximum number of `squeue` polls.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    task_ids = expand_array_spec(args.array_spec)
    (ROOT / "slurm" / "jobs").mkdir(parents=True, exist_ok=True)
    (ROOT / "slurm" / "logs").mkdir(parents=True, exist_ok=True)

    job_id, script_path = submit_array_job(
        job_name=args.job_name,
        partition=args.partition,
        array_spec=args.array_spec,
        sleep_seconds=args.sleep,
        cpus=args.cpus,
        mem=args.mem,
        time_limit=args.time_limit,
    )
    queue_snapshots = collect_squeue_snapshots(job_id, args.poll_interval, args.max_polls)

    task_summaries = []
    for task_id in task_ids:
        accounting = fetch_task_accounting(job_id, task_id)
        stdout_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}_{task_id}.out"
        stderr_path = ROOT / "slurm" / "logs" / f"{args.job_name}-{job_id}_{task_id}.err"
        task_summaries.append(
            {
                "task_id": task_id,
                "accounting": accounting,
                "stdout_path": str(stdout_path.relative_to(ROOT)),
                "stderr_path": str(stderr_path.relative_to(ROOT)),
                "stdout": read_text(stdout_path),
                "stderr": read_text(stderr_path),
            }
        )

    payload = {
        "job_id": job_id,
        "job_name": args.job_name,
        "partition": args.partition,
        "array_spec": args.array_spec,
        "task_count": len(task_ids),
        "task_ids": task_ids,
        "script_path": str(script_path.relative_to(ROOT)),
        "queue_snapshot_count": len(queue_snapshots),
        "queue_states_seen": sorted({row["State"] for snapshot in queue_snapshots for row in snapshot}),
        "queue_snapshots": queue_snapshots,
        "task_summaries": task_summaries,
    }
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
