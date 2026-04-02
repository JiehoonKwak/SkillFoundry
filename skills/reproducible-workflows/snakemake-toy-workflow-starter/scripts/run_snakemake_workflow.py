#!/usr/bin/env python3
"""Run the toy Snakemake workflow through the repo-managed Snakemake prefix."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = SKILL_ROOT / "examples"
SNAKEMAKE_PREFIX = ROOT / "slurm" / "envs" / "snakemake"


def snakemake_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{SNAKEMAKE_PREFIX / 'bin'}:{env['PATH']}"
    cache_home = ROOT / "scratch" / "snakemake-cache"
    cache_home.mkdir(parents=True, exist_ok=True)
    env["XDG_CACHE_HOME"] = str(cache_home)
    env["TMPDIR"] = str(ROOT / "scratch" / "tmp")
    Path(env["TMPDIR"]).mkdir(parents=True, exist_ok=True)
    return env


def run_command(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        check=True,
        capture_output=True,
        text=True,
        env=snakemake_env(),
        cwd=str(cwd or ROOT),
        timeout=240,
    )


def seed_workspace(workspace: Path) -> None:
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    for source in EXAMPLES.iterdir():
        destination = workspace / source.name
        if source.is_dir():
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument(
        "--workspace",
        type=Path,
        default=ROOT / "scratch" / "snakemake-toy-workflow",
        help="Workspace where the examples are copied and executed",
    )
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="Optional JSON output path for the run summary",
    )
    parser.add_argument(
        "--results-copy",
        type=Path,
        default=None,
        help="Optional directory to receive a copy of the results folder",
    )
    args = parser.parse_args()

    workspace = args.workspace.resolve()
    seed_workspace(workspace)

    version = run_command(["snakemake", "--version"])
    execution = run_command(
        [
            "snakemake",
            "--directory",
            str(workspace),
            "--snakefile",
            str(workspace / "Snakefile"),
            "--cores",
            "1",
            "--quiet",
            "all",
        ]
    )

    results_dir = workspace / "results"
    copied_path = results_dir / "copied.txt"
    summary_path = results_dir / "summary.json"
    summary_payload = json.loads(summary_path.read_text(encoding="utf-8"))
    copied_text = copied_path.read_text(encoding="utf-8").strip()

    if args.results_copy is not None:
        target = args.results_copy.resolve()
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(results_dir, target)

    payload = {
        "snakemake_version": version.stdout.strip(),
        "workspace": str(workspace),
        "results_dir": str(results_dir),
        "copied_text": copied_text,
        "result_files": sorted(path.name for path in results_dir.iterdir() if path.is_file()),
        "summary": summary_payload,
        "stdout": execution.stdout.strip().splitlines(),
        "stderr": execution.stderr.strip().splitlines(),
    }

    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
