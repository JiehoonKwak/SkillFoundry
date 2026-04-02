#!/usr/bin/env python3
"""Create a tiny repo and run local pre-commit hooks as a regression-testing starter."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PRECOMMIT_BIN = ROOT / "slurm" / "envs" / "maintenance" / "bin" / "pre-commit"

CONFIG_TEXT = """repos:
- repo: local
  hooks:
    - id: python-compile
      name: python-compile
      entry: python -m py_compile
      language: system
      types: [python]
    - id: json-validate
      name: json-validate
      entry: python -m json.tool
      language: system
      files: \\.json$
"""


def write_workspace(workspace: Path) -> dict:
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    (workspace / "src").mkdir()
    (workspace / "data").mkdir()
    (workspace / "src" / "toy_module.py").write_text(
        "def answer() -> int:\n    return 42\n",
        encoding="utf-8",
    )
    (workspace / "data" / "toy_metadata.json").write_text(
        json.dumps({"skill": "precommit-regression-testing-starter", "value": 42}, indent=2) + "\n",
        encoding="utf-8",
    )
    config_path = workspace / ".pre-commit-config.yaml"
    config_path.write_text(CONFIG_TEXT, encoding="utf-8")
    return {
        "config_path": str(config_path),
        "tracked_files": ["src/toy_module.py", "data/toy_metadata.json", ".pre-commit-config.yaml"],
    }


def run(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )


def run_precommit_regression(workspace: Path) -> dict:
    workspace_summary = write_workspace(workspace)
    run(["git", "init"], cwd=workspace)
    run(["git", "add", "."], cwd=workspace)
    env = os.environ.copy()
    env["PRE_COMMIT_HOME"] = str(workspace / ".pre-commit-cache")
    completed = run(
        [str(PRECOMMIT_BIN), "run", "--all-files", "--color", "never"],
        cwd=workspace,
        env=env,
    )
    stdout_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    return {
        "workspace": str(workspace),
        "config_path": workspace_summary["config_path"],
        "hook_ids": ["python-compile", "json-validate"],
        "tracked_file_count": len(workspace_summary["tracked_files"]),
        "returncode": int(completed.returncode),
        "stdout_tail": stdout_lines[-4:],
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True, help="Workspace path for the temporary repo.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if not PRECOMMIT_BIN.exists():
        raise SystemExit("pre-commit is not installed at slurm/envs/maintenance/bin/pre-commit")

    payload = run_precommit_regression(args.workspace)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
