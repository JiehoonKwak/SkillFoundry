#!/usr/bin/env python3
"""Validate experiment-only sc_skills on deterministic toy data."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "experiments" / "sc_skills" / "batch_design_manifest.json"
DEFAULT_SCRATCH_ROOT = ROOT / "scratch" / "sc_skills_validation"


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def select_tasks(manifest: dict, *, task_slugs: list[str], select_all: bool, limit: int | None) -> list[dict]:
    tasks = manifest["tasks"]
    if task_slugs:
        lookup = {task["task_slug"]: task for task in tasks}
        selected = [lookup[slug] for slug in task_slugs]
    elif select_all:
        selected = list(tasks)
    else:
        selected = list(tasks[:1])
    if limit is not None:
        selected = selected[:limit]
    return selected


def run_command(command: list[str], *, cwd: Path) -> dict:
    started = time.time()
    completed = subprocess.run(command, cwd=cwd, check=False, capture_output=True, text=True)
    return {
        "command": command,
        "cwd": str(cwd),
        "returncode": completed.returncode,
        "duration_seconds": round(time.time() - started, 3),
        "stdout_tail": completed.stdout.splitlines()[-20:],
        "stderr_tail": completed.stderr.splitlines()[-20:],
    }


def validate_task(task: dict, *, scratch_root: Path, run_local_tests: bool) -> dict:
    skill_root = ROOT / task["target_dir"]
    outdir = scratch_root / task["task_slug"]
    outdir.mkdir(parents=True, exist_ok=True)

    run_result = run_command(
        [sys.executable, str(skill_root / "scripts" / "run_exercise.py"), "--outdir", str(outdir)],
        cwd=ROOT,
    )
    validate_result = run_command(
        [sys.executable, str(skill_root / "scripts" / "validate_outputs.py"), "--outdir", str(outdir)],
        cwd=ROOT,
    )
    test_result = None
    if run_local_tests:
        test_result = run_command(
            [
                sys.executable,
                "-m",
                "unittest",
                "discover",
                "-s",
                str(skill_root / "tests"),
                "-p",
                "test_*.py",
            ],
            cwd=ROOT,
        )

    checks = [run_result["returncode"] == 0, validate_result["returncode"] == 0]
    if test_result is not None:
        checks.append(test_result["returncode"] == 0)

    missing_outputs: list[str] = []
    for relative_path in task["deliverables"]:
        if not (outdir / relative_path).exists():
            missing_outputs.append(relative_path)
    if not (outdir / "run_summary.json").exists():
        missing_outputs.append("run_summary.json")

    success = all(checks) and not missing_outputs
    payload = {
        "task_slug": task["task_slug"],
        "skill_root": str(skill_root),
        "outdir": str(outdir),
        "success": success,
        "missing_outputs": missing_outputs,
        "run_exercise": run_result,
        "validate_outputs": validate_result,
        "local_tests": test_result,
    }
    (outdir / "validation_result.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def render_markdown(results: list[dict], *, run_local_tests: bool) -> str:
    passed = sum(1 for item in results if item["success"])
    lines = [
        "# sc_skills Validation Report",
        "",
        "## Aggregate",
        "",
        "| Metric | Value |",
        "| --- | --- |",
        f"| Tasks validated | {len(results)} |",
        f"| Tasks passed | {passed} |",
        f"| Tasks failed | {len(results) - passed} |",
        f"| Local tests enabled | {run_local_tests} |",
        "",
        "## Per-task results",
        "",
        "| Task | Success | Run exercise | Validate outputs | Local tests | Missing outputs |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for item in results:
        test_rc = "-" if item["local_tests"] is None else str(item["local_tests"]["returncode"])
        missing = ", ".join(item["missing_outputs"]) if item["missing_outputs"] else "-"
        lines.append(
            f"| {item['task_slug']} | {item['success']} | {item['run_exercise']['returncode']} | {item['validate_outputs']['returncode']} | {test_rc} | {missing} |"
        )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", action="append", default=[], help="Task slug from the sc_skills manifest.")
    parser.add_argument("--all", action="store_true", help="Validate all experiment tasks.")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--scratch-root", type=Path, default=DEFAULT_SCRATCH_ROOT)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--markdown-out", type=Path, required=True)
    parser.add_argument("--run-local-tests", action="store_true", help="Also run unittest discovery inside each skill package.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    manifest = load_manifest()
    selected = select_tasks(manifest, task_slugs=args.task, select_all=args.all, limit=args.limit)
    args.scratch_root.mkdir(parents=True, exist_ok=True)

    results = [validate_task(task, scratch_root=args.scratch_root, run_local_tests=args.run_local_tests) for task in selected]
    payload = {
        "validated_tasks": [task["task_slug"] for task in selected],
        "total_tasks": len(results),
        "passed_tasks": sum(1 for item in results if item["success"]),
        "failed_tasks": [item["task_slug"] for item in results if not item["success"]],
        "run_local_tests": args.run_local_tests,
        "results": results,
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown_out.write_text(render_markdown(results, run_local_tests=args.run_local_tests), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if not payload["failed_tasks"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
