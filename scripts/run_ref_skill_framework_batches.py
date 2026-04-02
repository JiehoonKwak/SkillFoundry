#!/usr/bin/env python3
"""Run experiment sc-skill tasks through the framework in resumable batches."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "experiments" / "sc_skills" / "batch_design_manifest.json"
DEFAULT_STATE_ROOT = ROOT / "scratch" / "framework" / "sc-skill-batches"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def ordered_task_slugs(manifest: dict) -> list[str]:
    return [task["task_slug"] for task in manifest["tasks"]]


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[index:index + size] for index in range(0, len(items), size)]


def state_paths(state_root: Path, label: str) -> dict[str, Path]:
    base = state_root / label
    base.mkdir(parents=True, exist_ok=True)
    return {
        "base": base,
        "state": base / "state.json",
        "status": base / "status.json",
        "summary": base / "summary.md",
        "log": base / "batch.log",
        "tasks": base / "tasks",
    }


def logger_for(label: str, log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"sc-skill-batches.{label}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def init_state(
    *,
    label: str,
    batch_size: int,
    verification_mode: str,
    framework_model: str | None,
    framework_reasoning_effort: str | None,
    framework_profile: str | None,
    framework_codex_bin: str | None,
    framework_full_auto: bool,
    framework_extra_context: str | None,
) -> dict:
    manifest = load_manifest()
    return {
        "label": label,
        "started_at": utc_now(),
        "updated_at": utc_now(),
        "status": "initialized",
        "batch_size": batch_size,
        "verification_mode": verification_mode,
        "framework_model": framework_model,
        "framework_reasoning_effort": framework_reasoning_effort,
        "framework_profile": framework_profile,
        "framework_codex_bin": framework_codex_bin,
        "framework_full_auto": framework_full_auto,
        "framework_extra_context": framework_extra_context,
        "task_order": ordered_task_slugs(manifest),
        "completed_tasks": [],
        "failed_tasks": [],
        "task_results": {},
        "batches": [],
    }


def load_or_init_state(
    *,
    label: str,
    batch_size: int,
    verification_mode: str,
    framework_model: str | None,
    framework_reasoning_effort: str | None,
    framework_profile: str | None,
    framework_codex_bin: str | None,
    framework_full_auto: bool,
    framework_extra_context: str | None,
    paths: dict[str, Path],
) -> dict:
    if paths["state"].exists():
        payload = json.loads(paths["state"].read_text(encoding="utf-8"))
        payload["updated_at"] = utc_now()
        payload["framework_model"] = framework_model
        payload["framework_reasoning_effort"] = framework_reasoning_effort
        payload["framework_profile"] = framework_profile
        payload["framework_codex_bin"] = framework_codex_bin
        payload["framework_full_auto"] = framework_full_auto
        payload["framework_extra_context"] = framework_extra_context
        return payload
    return init_state(
        label=label,
        batch_size=batch_size,
        verification_mode=verification_mode,
        framework_model=framework_model,
        framework_reasoning_effort=framework_reasoning_effort,
        framework_profile=framework_profile,
        framework_codex_bin=framework_codex_bin,
        framework_full_auto=framework_full_auto,
        framework_extra_context=framework_extra_context,
    )


def pending_tasks(state: dict) -> list[str]:
    completed = set(state.get("completed_tasks", []))
    return [slug for slug in state["task_order"] if slug not in completed]


def status_payload(state: dict) -> dict:
    pending = pending_tasks(state)
    return {
        "label": state["label"],
        "status": state["status"],
        "started_at": state.get("started_at"),
        "updated_at": utc_now(),
        "batch_size": state["batch_size"],
        "verification_mode": state["verification_mode"],
        "framework_model": state.get("framework_model"),
        "framework_reasoning_effort": state.get("framework_reasoning_effort"),
        "total_tasks": len(state["task_order"]),
        "completed_tasks": len(state.get("completed_tasks", [])),
        "failed_tasks": len(state.get("failed_tasks", [])),
        "pending_tasks": len(pending),
        "pending_preview": pending[:10],
        "batches_completed": sum(1 for batch in state.get("batches", []) if batch.get("status") == "completed"),
        "batches_started": len(state.get("batches", [])),
        "last_batch_index": state.get("batches", [{}])[-1].get("batch_index") if state.get("batches") else None,
    }


def write_summary_markdown(path: Path, state: dict, status: dict) -> None:
    lines = [
        f"# sc_skill Batch Summary: {state['label']}",
        "",
        "## Aggregate",
        "",
        "| Metric | Value |",
        "| --- | --- |",
        f"| Status | {status['status']} |",
        f"| Total tasks | {status['total_tasks']} |",
        f"| Completed tasks | {status['completed_tasks']} |",
        f"| Failed tasks | {status['failed_tasks']} |",
        f"| Pending tasks | {status['pending_tasks']} |",
        f"| Batches started | {status['batches_started']} |",
        f"| Batches completed | {status['batches_completed']} |",
        f"| Verification mode | {status['verification_mode']} |",
        f"| Framework model | {status.get('framework_model') or 'default'} |",
        f"| Framework reasoning effort | {status.get('framework_reasoning_effort') or 'default'} |",
        "",
        "## Recent batches",
        "",
        "| Batch | Status | Tasks |",
        "| --- | --- | --- |",
    ]
    for batch in state.get("batches", [])[-10:]:
        lines.append(f"| {batch['batch_index']} | {batch['status']} | {', '.join(batch.get('task_slugs', []))} |")
    lines.extend(["", "## Pending preview", ""])
    for slug in status["pending_preview"]:
        lines.append(f"- {slug}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_state(paths: dict[str, Path], state: dict) -> dict:
    state["updated_at"] = utc_now()
    paths["state"].write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status = status_payload(state)
    paths["status"].write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_summary_markdown(paths["summary"], state, status)
    return status


def framework_command(
    task_slug: str,
    verification_mode: str,
    *,
    framework_model: str | None,
    framework_reasoning_effort: str | None,
    framework_profile: str | None,
    framework_codex_bin: str | None,
    framework_full_auto: bool,
    framework_extra_context: str | None,
) -> list[str]:
    command = [
        sys.executable,
        "scripts/batch_design_experiment_skills.py",
        "--task",
        task_slug,
        "--run-framework",
        "--verification-mode",
        verification_mode,
    ]
    if framework_model:
        command.extend(["--framework-model", framework_model])
    if framework_reasoning_effort:
        command.extend(["--framework-reasoning-effort", framework_reasoning_effort])
    if framework_profile:
        command.extend(["--framework-profile", framework_profile])
    if framework_codex_bin:
        command.extend(["--framework-codex-bin", framework_codex_bin])
    if framework_full_auto:
        command.append("--framework-full-auto")
    if framework_extra_context:
        command.extend(["--framework-extra-context", framework_extra_context])
    return command


def validation_command(task_slug: str, *, task_dir: Path) -> list[str]:
    return [
        sys.executable,
        "scripts/validate_sc_skill_experiments.py",
        "--task",
        task_slug,
        "--scratch-root",
        str(task_dir / "artifacts"),
        "--json-out",
        str(task_dir / "validation_report.json"),
        "--markdown-out",
        str(task_dir / "validation_report.md"),
        "--run-local-tests",
    ]


def run_task(
    task_slug: str,
    *,
    verification_mode: str,
    framework_model: str | None,
    framework_reasoning_effort: str | None,
    framework_profile: str | None,
    framework_codex_bin: str | None,
    framework_full_auto: bool,
    framework_extra_context: str | None,
    task_log_root: Path,
) -> dict:
    task_log_root.mkdir(parents=True, exist_ok=True)
    task_dir = task_log_root / task_slug
    task_dir.mkdir(parents=True, exist_ok=True)
    console_path = task_dir / "console.log"
    framework_stdout_path = task_dir / "framework.stdout.log"
    framework_stderr_path = task_dir / "framework.stderr.log"
    validation_stdout_path = task_dir / "validation.stdout.log"
    validation_stderr_path = task_dir / "validation.stderr.log"
    started = time.time()
    framework_cmd = framework_command(
        task_slug,
        verification_mode,
        framework_model=framework_model,
        framework_reasoning_effort=framework_reasoning_effort,
        framework_profile=framework_profile,
        framework_codex_bin=framework_codex_bin,
        framework_full_auto=framework_full_auto,
        framework_extra_context=framework_extra_context,
    )
    framework_completed = subprocess.run(
        framework_cmd,
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    framework_stdout_path.write_text(framework_completed.stdout, encoding="utf-8")
    framework_stderr_path.write_text(framework_completed.stderr, encoding="utf-8")
    validation_completed = None
    if framework_completed.returncode == 0:
        validation_completed = subprocess.run(
            validation_command(task_slug, task_dir=task_dir),
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        validation_stdout_path.write_text(validation_completed.stdout, encoding="utf-8")
        validation_stderr_path.write_text(validation_completed.stderr, encoding="utf-8")
    console_path.write_text(
        "\n".join(
            [
                "# Framework command",
                " ".join(framework_cmd),
                "",
                "# Framework STDOUT",
                framework_completed.stdout.rstrip(),
                "",
                "# Framework STDERR",
                framework_completed.stderr.rstrip(),
                "",
                "# Validation command",
                " ".join(validation_command(task_slug, task_dir=task_dir)) if validation_completed is not None else "skipped",
                "",
                "# Validation STDOUT",
                "" if validation_completed is None else validation_completed.stdout.rstrip(),
                "",
                "# Validation STDERR",
                "" if validation_completed is None else validation_completed.stderr.rstrip(),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    framework_stdout_lines = framework_completed.stdout.splitlines()
    framework_stderr_lines = framework_completed.stderr.splitlines()
    validation_stdout_lines = [] if validation_completed is None else validation_completed.stdout.splitlines()
    validation_stderr_lines = [] if validation_completed is None else validation_completed.stderr.splitlines()
    aggregate_returncode = framework_completed.returncode if validation_completed is None else max(framework_completed.returncode, validation_completed.returncode)
    duration = round(time.time() - started, 3)
    payload = {
        "task_slug": task_slug,
        "returncode": aggregate_returncode,
        "duration_seconds": duration,
        "started_at": utc_now(),
        "framework_returncode": framework_completed.returncode,
        "framework_stdout_tail": framework_stdout_lines[-20:],
        "framework_stderr_tail": framework_stderr_lines[-20:],
        "validation_returncode": None if validation_completed is None else validation_completed.returncode,
        "validation_stdout_tail": validation_stdout_lines[-20:],
        "validation_stderr_tail": validation_stderr_lines[-20:],
        "task_dir": str(task_dir),
        "command": framework_cmd,
        "validation_command": None if validation_completed is None else validation_command(task_slug, task_dir=task_dir),
    }
    if framework_stdout_lines:
        try:
            payload["result"] = json.loads(framework_completed.stdout)
        except json.JSONDecodeError:
            payload["result"] = None
    else:
        payload["result"] = None
    (task_dir / "result.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="Run sc-skill framework batches.")
    run.add_argument("--label", default="sc-skill-framework-batches")
    run.add_argument("--batch-size", type=int, default=4)
    run.add_argument("--verification-mode", default="none", choices=["none", "validate", "standard", "full", "audit"])
    run.add_argument("--state-root", type=Path, default=DEFAULT_STATE_ROOT)
    run.add_argument("--max-batches", type=int, default=None, help="Optional number of batches to run this invocation.")
    run.add_argument("--max-runtime-minutes", type=int, default=450)
    run.add_argument("--stop-buffer-minutes", type=int, default=20)
    run.add_argument("--print-plan", action="store_true", help="Print the current batch plan and exit.")
    run.add_argument("--framework-model", default=None)
    run.add_argument("--framework-reasoning-effort", default=None)
    run.add_argument("--framework-profile", default=None)
    run.add_argument("--framework-codex-bin", default=None)
    run.add_argument("--framework-full-auto", action="store_true")
    run.add_argument("--framework-extra-context", default=None)

    status = subparsers.add_parser("status", help="Show the latest batch status.")
    status.add_argument("--label", required=True)
    status.add_argument("--state-root", type=Path, default=DEFAULT_STATE_ROOT)

    return parser


def cmd_status(args: argparse.Namespace) -> int:
    paths = state_paths(args.state_root, args.label)
    if not paths["status"].exists():
        print(json.dumps({"label": args.label, "status": "missing"}, indent=2))
        return 1
    print(paths["status"].read_text(encoding="utf-8"))
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    paths = state_paths(args.state_root, args.label)
    log = logger_for(args.label, paths["log"])
    state = load_or_init_state(
        label=args.label,
        batch_size=args.batch_size,
        verification_mode=args.verification_mode,
        framework_model=args.framework_model,
        framework_reasoning_effort=args.framework_reasoning_effort,
        framework_profile=args.framework_profile,
        framework_codex_bin=args.framework_codex_bin,
        framework_full_auto=args.framework_full_auto,
        framework_extra_context=args.framework_extra_context,
        paths=paths,
    )
    remaining = pending_tasks(state)
    batch_plan = chunked(remaining, state["batch_size"])

    if args.print_plan:
        print(
            json.dumps(
                {
                    "label": args.label,
                    "batch_size": state["batch_size"],
                    "remaining_tasks": remaining,
                    "batch_plan": batch_plan,
                },
                indent=2,
            )
        )
        return 0

    log.info("Starting sc-skill batch runner: label=%s batch_size=%s pending=%s", args.label, state["batch_size"], len(remaining))
    started_at = time.time()
    max_runtime_seconds = max(0, args.max_runtime_minutes * 60)
    stop_buffer_seconds = max(0, args.stop_buffer_minutes * 60)
    batches_run = 0

    for batch_index, task_slugs in enumerate(batch_plan, start=1):
        elapsed = time.time() - started_at
        if max_runtime_seconds and elapsed >= max_runtime_seconds - stop_buffer_seconds:
            state["status"] = "checkpointed"
            write_state(paths, state)
            log.info("Stopping before batch %s to respect runtime budget.", batch_index)
            return 0
        if args.max_batches is not None and batches_run >= args.max_batches:
            state["status"] = "checkpointed"
            write_state(paths, state)
            log.info("Stopping after requested max-batches=%s.", args.max_batches)
            return 0

        batch_record = {"batch_index": batch_index, "task_slugs": task_slugs, "status": "running", "results": []}
        state["batches"].append(batch_record)
        state["status"] = "running"
        write_state(paths, state)
        log.info("Starting batch %s with tasks: %s", batch_index, ", ".join(task_slugs))

        for task_slug in task_slugs:
            log.info("Starting task %s (verification=%s model=%s effort=%s)", task_slug, state["verification_mode"], state.get("framework_model") or "default", state.get("framework_reasoning_effort") or "default")
            result = run_task(
                task_slug,
                verification_mode=args.verification_mode,
                framework_model=state.get("framework_model"),
                framework_reasoning_effort=state.get("framework_reasoning_effort"),
                framework_profile=state.get("framework_profile"),
                framework_codex_bin=state.get("framework_codex_bin"),
                framework_full_auto=bool(state.get("framework_full_auto")),
                framework_extra_context=state.get("framework_extra_context"),
                task_log_root=paths["tasks"],
            )
            batch_record["results"].append(result)
            state["task_results"][task_slug] = result
            if result["returncode"] == 0:
                state["completed_tasks"].append(task_slug)
                log.info("Completed task %s successfully in %.3fs", task_slug, result["duration_seconds"])
            else:
                state["failed_tasks"].append(task_slug)
                log.warning("Task %s failed in %.3fs (framework=%s validation=%s)", task_slug, result["duration_seconds"], result["framework_returncode"], result["validation_returncode"])
            write_state(paths, state)

        batch_record["status"] = "completed"
        batches_run += 1
        log.info("Completed batch %s", batch_index)
        write_state(paths, state)

    state["status"] = "completed"
    write_state(paths, state)
    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "status":
        return cmd_status(args)
    return cmd_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
