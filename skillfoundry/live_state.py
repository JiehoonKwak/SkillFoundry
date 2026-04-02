"""Live run state helpers for framework dashboard polling."""

from __future__ import annotations

import os
import sys
import threading
from datetime import datetime
from pathlib import Path

from .repository import write_json


class FrameworkRunStateWriter:
    """Persist a lightweight, poll-friendly snapshot of the currently active run."""

    COLORS = {
        "reset": "\033[0m",
        "dim": "\033[2m",
        "info": "\033[36m",
        "command": "\033[35m",
        "stdout": "\033[37m",
        "stderr": "\033[2m",
        "success": "\033[32m",
        "warning": "\033[33m",
        "error": "\033[31m",
        "stage": "\033[34m",
    }

    def __init__(self, state_dir: str | Path, *, log_limit: int = 80, emit_console: bool | None = None) -> None:
        self.state_dir = Path(state_dir)
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.path = self.state_dir / "latest_active_run.json"
        self.log_limit = log_limit
        if emit_console is None:
            emit_console = os.environ.get("SCISKILL_FRAMEWORK_EMIT_CONSOLE", "1") != "0"
        self.emit_console = emit_console and os.environ.get("SCISKILL_FRAMEWORK_SILENT") != "1"
        self.use_color = self.emit_console and (
            sys.stderr.isatty() or os.environ.get("FORCE_COLOR") == "1"
        )
        self._lock = threading.RLock()
        self.state: dict = {}

    def _now(self) -> str:
        return datetime.now().isoformat(timespec="seconds")

    def _flush(self) -> None:
        with self._lock:
            write_json(self.path, self.state)

    def _console(self, prefix: str, text: str, *, color: str | None = None) -> None:
        if not self.emit_console:
            return
        if self.use_color and color in self.COLORS:
            message = f"{self.COLORS[color]}{prefix}{self.COLORS['reset']} {text}"
        else:
            message = f"{prefix} {text}"
        print(message, file=sys.stderr, flush=True)

    def start_run(
        self,
        *,
        mode: str,
        run_label: str,
        run_dir: str,
        verification_mode: str,
        steps: list[dict],
        selected_skills: list[str] | None = None,
    ) -> None:
        self.state = {
            "active": True,
            "status": "running",
            "mode": mode,
            "run_label": run_label,
            "run_dir": run_dir,
            "verification_mode": verification_mode,
            "started_at": self._now(),
            "updated_at": self._now(),
            "ended_at": None,
            "current_step": None,
            "selected_skills": selected_skills or [],
            "steps": steps,
            "log_tail": [],
        }
        self._flush()
        selected_count = len(selected_skills or [])
        self._console(
            "[RUN]",
            f"{mode} started: {run_label} | skills={selected_count} | verification={verification_mode} | run_dir={run_dir}",
            color="stage",
        )

    def append_log(self, *, kind: str, label: str, text: str) -> None:
        if not text:
            return
        with self._lock:
            logs = self.state.setdefault("log_tail", [])
            logs.append(
                {
                    "at": self._now(),
                    "kind": kind,
                    "label": label,
                    "text": text,
                }
            )
            if len(logs) > self.log_limit:
                self.state["log_tail"] = logs[-self.log_limit :]
            self.state["updated_at"] = self._now()
            self._flush()
        prefix = {
            "info": "[INFO]",
            "command": "[CMD ]",
            "stdout": "[OUT ]",
            "stderr": "[STDERR]",
            "success": "[ OK ]",
            "warning": "[WARN]",
            "error": "[ERR ]",
        }.get(kind, "[LOG ]")
        self._console(prefix, f"{label}: {text}", color=kind if kind in self.COLORS else None)

    def append_lines(self, *, kind: str, label: str, lines: list[str], limit: int = 8) -> None:
        trimmed = [line for line in (lines or []) if line][-limit:]
        for line in trimmed:
            self.append_log(kind=kind, label=label, text=line)

    def mark_step_running(
        self,
        step_name: str,
        *,
        kind: str,
        summary: str | None = None,
        detail: dict | None = None,
    ) -> None:
        with self._lock:
            for step in self.state.get("steps", []):
                if step["name"] == step_name:
                    step["status"] = "running"
                    if summary:
                        step["summary"] = summary
                    if detail:
                        step["detail"] = detail
                    break
            self.state["current_step"] = {
                "name": step_name,
                "kind": kind,
                "summary": summary,
                **(detail or {}),
            }
            self.state["updated_at"] = self._now()
            self._flush()
        self._console(
            "[STEP]",
            f"{step_name}: {summary or 'running'}",
            color="stage",
        )

    def mark_step_complete(
        self,
        step_name: str,
        *,
        summary: str,
        health: str = "clean",
        blockers: list[str] | None = None,
    ) -> None:
        blockers = blockers or []
        status = "attention" if blockers or health != "clean" else "completed"
        with self._lock:
            for step in self.state.get("steps", []):
                if step["name"] == step_name:
                    step["status"] = status
                    step["health"] = health
                    step["summary"] = summary
                    step["blockers"] = blockers
                    break
            self.state["current_step"] = {
                "name": step_name,
                "kind": "completed",
                "summary": summary,
            }
            self.state["updated_at"] = self._now()
            self._flush()
        self._console(
            "[DONE]" if status == "completed" else "[WARN]",
            f"{step_name}: {summary}",
            color="success" if status == "completed" else "warning",
        )

    def finish(self, *, status: str, summary: str) -> None:
        with self._lock:
            self.state["active"] = False
            self.state["status"] = status
            self.state["ended_at"] = self._now()
            self.state["updated_at"] = self._now()
            self.state["current_step"] = {
                "name": "completed",
                "kind": "completed",
                "summary": summary,
            }
        self.append_log(
            kind="success" if status == "completed" else "warning",
            label="run",
            text=summary,
        )
        self._flush()
