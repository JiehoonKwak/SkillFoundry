#!/usr/bin/env python3
"""Run the local WDL greeting example with miniwdl+udocker and summarize outputs."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
HOME_UDOCKER = Path.home() / ".udocker"
UDOCKER_CONTAINER_NAME = "miniwdlubuntu"
UDOCKER_IMAGE = "ubuntu:20.04"
MINIWDL = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "miniwdl"
UDOCKER = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "udocker"
EXAMPLE_WDL = SKILL_DIR / "examples" / "hello.wdl"


def build_inputs_file(workspace: Path, name: str) -> Path:
    inputs = {"name": name}
    path = workspace / "hello-inputs.json"
    path.write_text(json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def build_config_file(workspace: Path) -> Path:
    path = workspace / "miniwdl-local.cfg"
    path.write_text(
        "\n".join(
            [
                "[scheduler]",
                "container_backend = udocker",
                "",
                "[task_runtime]",
                "defaults = {",
                f'        "docker": "{UDOCKER_CONTAINER_NAME}"',
                "    }",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def seed_udocker_workspace(workspace: Path) -> Path:
    path = workspace / ".udocker"
    if path.exists():
        return path
    if not HOME_UDOCKER.exists():
        raise SystemExit(f"Missing udocker cache directory: {HOME_UDOCKER}")
    shutil.copytree(HOME_UDOCKER, path, symlinks=True)
    return path


def run_udocker(args: list[str], *, udocker_dir: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    run_env = dict(os.environ if env is None else env)
    run_env["UDOCKER_DIR"] = str(udocker_dir)
    return subprocess.run(
        [str(UDOCKER), "--allow-root", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
        env=run_env,
    )


def ensure_udocker_container(udocker_dir: Path) -> None:
    listing = run_udocker(["ps", "-m", "-s"], udocker_dir=udocker_dir)
    if UDOCKER_CONTAINER_NAME in listing.stdout:
        return
    create_env = dict(os.environ)
    create_env["UDOCKER_DEFAULT_EXECUTION_MODE"] = "F3"
    run_udocker(
        ["create", f"--name={UDOCKER_CONTAINER_NAME}", UDOCKER_IMAGE],
        udocker_dir=udocker_dir,
        env=create_env,
    )
    run_udocker(["setup", "--execmode=F4", UDOCKER_CONTAINER_NAME], udocker_dir=udocker_dir)


def build_udocker_shim(workspace: Path) -> Path:
    path = workspace / "udocker-miniwdl.sh"
    path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f'real_udocker="{UDOCKER}"',
                'for arg in "$@"; do',
                '  case "$arg" in',
                '    --version|-V|version)',
                '      exec "$real_udocker" "$@"',
                '      ;;',
                '  esac',
                'done',
                'cmd=""',
                'for arg in "$@"; do',
                '  if [[ "$arg" != -* ]]; then',
                '    cmd="$arg"',
                '    break',
                '  fi',
                'done',
                'if [[ "$cmd" == "pull" ]]; then',
                '  exit 0',
                'fi',
                'exec "$real_udocker" "$@"',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def run_command(name: str, workspace: Path) -> dict:
    if not MINIWDL.exists():
        raise SystemExit(f"Missing miniwdl executable: {MINIWDL}")
    if not UDOCKER.exists():
        raise SystemExit(f"Missing udocker executable: {UDOCKER}")
    workspace = workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    udocker_dir = seed_udocker_workspace(workspace)
    ensure_udocker_container(udocker_dir)
    summary_json = (workspace / "miniwdl-outputs.json").resolve()
    inputs_path = build_inputs_file(workspace, name).resolve()
    config_path = build_config_file(workspace).resolve()
    udocker_shim = build_udocker_shim(workspace).resolve()
    env = dict(os.environ)
    env["MINIWDL__SCHEDULER__CONTAINER_BACKEND"] = "udocker"
    env["UDOCKER_DIR"] = str(udocker_dir)
    env["MINIWDL__UDOCKER__EXE"] = json.dumps([str(udocker_shim), "--allow-root"])
    completed = subprocess.run(
        [
            str(MINIWDL),
            "run",
            str(EXAMPLE_WDL),
            "--input",
            str(inputs_path),
            "--dir",
            str(workspace),
            "--cfg",
            str(config_path),
            "-o",
            str(summary_json),
            "--error-json",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
        env=env,
    )
    outputs = json.loads(summary_json.read_text(encoding="utf-8"))
    output_file = Path(outputs["outputs"]["hello_workflow.greeting_file"])
    return {
        "name": name,
        "miniwdl_version": subprocess.run(
            [str(MINIWDL), "--version"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout.strip().split()[-1],
        "wdl_document": str(EXAMPLE_WDL.relative_to(ROOT)),
        "config_file": str(config_path),
        "miniwdl_run_dir": outputs["dir"],
        "output_file": str(output_file),
        "output_text": outputs["outputs"]["hello_workflow.greeting_text"].strip(),
        "greeting_text": outputs["outputs"]["hello_workflow.greeting_text"].strip(),
        "raw_stdout": completed.stdout.strip(),
        "outputs": outputs,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", default="WDL", help="Name to render in the greeting")
    parser.add_argument("--workspace", type=Path, default=None, help="Optional workspace directory")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output")
    args = parser.parse_args()

    if args.workspace is None:
        with tempfile.TemporaryDirectory(prefix="wdl-hello-") as tmp_dir:
            payload = run_command(args.name, Path(tmp_dir))
    else:
        payload = run_command(args.name, args.workspace)

    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
