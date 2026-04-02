#!/usr/bin/env python3
"""Run a minimal Nextflow hello workflow through the repo-local toolchain."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_PREFIX = ROOT / "slurm" / "envs" / "nextflow-tools"


def nextflow_env() -> dict[str, str]:
    env = os.environ.copy()
    env["JAVA_HOME"] = str(NEXTFLOW_PREFIX)
    env["PATH"] = f"{NEXTFLOW_PREFIX / 'bin'}:{env['PATH']}"
    env.setdefault("NXF_ANSI_LOG", "false")
    return env


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        check=True,
        capture_output=True,
        text=True,
        env=nextflow_env(),
        cwd=ROOT,
        timeout=600,
    )


def write_slurm_config(path: Path, partition: str, cpus: int, memory: str, time_limit: str) -> None:
    path.write_text(
        "\n".join(
            [
                "process.executor = 'slurm'",
                f"process.queue = '{partition}'",
                f"process.cpus = {cpus}",
                f"process.memory = '{memory}'",
                f"process.time = '{time_limit}'",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def parse_trace(trace_path: Path) -> list[dict[str, str]]:
    if not trace_path.exists():
        return []
    with trace_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def collect_output_files(out_dir: Path, attempts: int = 10, delay: float = 1.0) -> list[dict[str, str]]:
    for attempt in range(attempts):
        files = []
        for path in sorted(out_dir.glob("*.txt")):
            files.append({"name": path.name, "text": path.read_text(encoding="utf-8").strip()})
        if files or attempt + 1 == attempts:
            return files
        time.sleep(delay)
    return []


def collect_trace_rows(trace_path: Path | None, attempts: int = 10, delay: float = 1.0) -> list[dict[str, str]]:
    if trace_path is None:
        return []
    for attempt in range(attempts):
        rows = parse_trace(trace_path)
        if rows or attempt + 1 == attempts:
            return rows
        time.sleep(delay)
    return []


def parse_sacct_table(text: str) -> dict[str, str]:
    lines = [line for line in text.strip().splitlines() if line.strip()]
    if len(lines) < 2:
        return {}
    headers = lines[0].split("|")
    values = lines[1].split("|")
    return dict(zip(headers, values, strict=True))


def fetch_accounting(job_id: str, attempts: int = 6) -> dict[str, str]:
    transitional_states = {"", "PENDING", "RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED"}
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
        state = parsed.get("State", "") if parsed else ""
        if parsed and state not in transitional_states:
            return parsed
        if attempt + 1 < attempts:
            time.sleep(1)
    raise RuntimeError(f"sacct did not return accounting for job {job_id}")


def is_tmp_path(path: Path) -> bool:
    try:
        return path.is_relative_to(Path("/tmp"))
    except AttributeError:
        return str(path).startswith("/tmp/")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument(
        "--pipeline",
        type=Path,
        default=SKILL_ROOT / "examples" / "main.nf",
        help="Path to the local Nextflow script",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=SKILL_ROOT / "assets" / "latest-results",
        help="Directory where published outputs will be written",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=ROOT / "scratch" / "nextflow-hello-work",
        help="Directory for Nextflow work files",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional JSON output path for the run summary",
    )
    parser.add_argument(
        "--executor",
        choices=("local", "slurm"),
        default="local",
        help="Execution backend for the smoke workflow",
    )
    parser.add_argument("--partition", default="cpu", help="Slurm partition for cluster-backed runs")
    parser.add_argument("--cpus", type=int, default=1, help="Requested CPUs for Slurm-backed runs")
    parser.add_argument("--mem", default="512 MB", help="Requested memory for Slurm-backed runs")
    parser.add_argument("--time-limit", default="5 min", help="Requested time limit for Slurm-backed runs")
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    work_dir = args.work_dir.resolve()
    if args.executor == "slurm" and (is_tmp_path(out_dir) or is_tmp_path(work_dir)):
        raise SystemExit("Slurm-backed Nextflow runs require shared filesystem paths such as scratch/ or repo-local directories, not /tmp.")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    info = run_command(["nextflow", "info"])
    command = [
        "nextflow",
        "run",
        str(args.pipeline.resolve()),
        "-ansi-log",
        "false",
        "-work-dir",
        str(work_dir),
        "--outdir",
        str(out_dir),
    ]
    trace_path = None
    if args.executor == "slurm":
        config_path = work_dir / "slurm-executor.config"
        trace_path = work_dir / "trace.txt"
        write_slurm_config(config_path, args.partition, args.cpus, args.mem, args.time_limit)
        command.extend(["-c", str(config_path), "-with-trace", str(trace_path)])
    run = run_command(command)

    if args.executor == "slurm":
        files = collect_output_files(out_dir, attempts=60, delay=2.0)
        trace_rows = collect_trace_rows(trace_path, attempts=60, delay=2.0)
    else:
        files = collect_output_files(out_dir, attempts=5)
        trace_rows = collect_trace_rows(trace_path, attempts=5)
    slurm_jobs = []
    if args.executor == "slurm":
        for row in trace_rows:
            native_id = row.get("native_id")
            if native_id:
                slurm_jobs.append(fetch_accounting(native_id))

    payload = {
        "executor": args.executor,
        "pipeline": str(args.pipeline.resolve()),
        "out_dir": str(out_dir),
        "work_dir": str(work_dir),
        "nextflow_info": info.stdout.strip().splitlines(),
        "stdout_tail": run.stdout.strip().splitlines()[-12:],
        "files": files,
        "trace_path": str(trace_path) if trace_path is not None else None,
        "trace_rows": trace_rows,
        "slurm_jobs": slurm_jobs,
    }
    if args.executor == "slurm":
        payload["partition"] = args.partition
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
