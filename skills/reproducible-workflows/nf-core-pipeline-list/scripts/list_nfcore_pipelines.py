#!/usr/bin/env python3
"""List nf-core pipelines through the repo-local environment and clean broken JSON."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
NEXTFLOW_PREFIX = ROOT / "slurm" / "envs" / "nextflow-tools"
ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")


def nfcore_env() -> dict[str, str]:
    env = os.environ.copy()
    env["JAVA_HOME"] = str(NEXTFLOW_PREFIX)
    env["PATH"] = f"{NEXTFLOW_PREFIX / 'bin'}:{env['PATH']}"
    env["NO_COLOR"] = "1"
    env["TERM"] = "dumb"
    return env


def repair_json(text: str) -> str:
    text = ANSI_RE.sub("", text).replace("\r", "")
    repaired: list[str] = []
    in_string = False
    escaped = False
    for char in text:
        if in_string:
            if escaped:
                repaired.append(char)
                escaped = False
                continue
            if char == "\\":
                repaired.append(char)
                escaped = True
                continue
            if char == "\"":
                repaired.append(char)
                in_string = False
                continue
            if char in "\n\t":
                repaired.append(" ")
                continue
            repaired.append(char)
            continue
        if char == "\"":
            in_string = True
        repaired.append(char)
    return "".join(repaired)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sort", default="pulled", choices=["release", "pulled", "name", "stars"])
    parser.add_argument("--limit", type=int, default=5, help="Maximum number of remote workflows to return")
    parser.add_argument(
        "--fixture",
        type=Path,
        default=None,
        help="Optional captured nf-core stdout fixture to parse instead of running the CLI",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.fixture is None:
        completed = subprocess.run(
            ["nf-core", "pipelines", "list", "--sort", args.sort, "--json"],
            cwd=ROOT,
            env=nfcore_env(),
            capture_output=True,
            text=True,
            check=True,
            timeout=180,
        )
        raw_text = completed.stdout
    else:
        raw_text = args.fixture.read_text(encoding="utf-8")

    payload = json.loads(repair_json(raw_text))
    local = [
        {
            "full_name": item["full_name"],
            "local_path": item["local_path"],
            "last_pull_date": item.get("last_pull_date"),
        }
        for item in payload.get("local_workflows", [])
    ]
    remote = []
    for item in payload.get("remote_workflows", [])[: args.limit]:
        latest_release = item.get("releases", [{}])[0]
        remote.append(
            {
                "name": item["name"],
                "full_name": item["full_name"],
                "archived": item["archived"],
                "stargazers_count": item["stargazers_count"],
                "latest_release": latest_release.get("tag_name"),
                "latest_release_date": latest_release.get("published_at"),
                "nextflow_version": latest_release.get("nextflow_version"),
                "nf_core_version": latest_release.get("nf_core_version"),
            }
        )

    summary = {
        "sort": args.sort,
        "counts": {
            "local_workflows": len(payload.get("local_workflows", [])),
            "remote_workflows": len(payload.get("remote_workflows", [])),
        },
        "local_workflows": local,
        "remote_workflows": remote,
    }
    text = json.dumps(summary, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
