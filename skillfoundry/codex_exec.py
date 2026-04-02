"""`codex exec` wrapper used by the framework."""

from __future__ import annotations

import json
import subprocess
import threading
import time
from pathlib import Path

from .models import CodexInvocationResult
from .repository import write_json


class CodexExecRunner:
    """Run structured `codex exec` stages and capture artifacts."""

    def __init__(
        self,
        codex_bin: str = "codex",
        model: str | None = None,
        reasoning_effort: str | None = None,
        profile: str | None = None,
        full_auto: bool = False,
        config_overrides: list[str] | None = None,
    ) -> None:
        self.codex_bin = codex_bin
        self.model = model
        self.reasoning_effort = reasoning_effort
        self.profile = profile
        self.full_auto = full_auto
        self.config_overrides = list(config_overrides or [])

    def run(
        self,
        *,
        prompt: str,
        output_schema: dict,
        working_dir: Path,
        artifact_dir: Path,
        model: str | None = None,
        reasoning_effort: str | None = None,
        profile: str | None = None,
        config_overrides: list[str] | None = None,
        timeout: int = 3600,
        progress_callback=None,
        heartbeat_seconds: int = 20,
        max_retries: int = 2,
    ) -> CodexInvocationResult:
        artifact_dir.mkdir(parents=True, exist_ok=True)
        prompt_path = artifact_dir / "prompt.md"
        schema_path = artifact_dir / "schema.json"
        stdout_path = artifact_dir / "stdout.txt"
        stderr_path = artifact_dir / "stderr.txt"
        last_message_path = artifact_dir / "last_message.json"
        attempts_path = artifact_dir / "attempts.json"

        prompt_path.write_text(prompt, encoding="utf-8")
        write_json(schema_path, output_schema)

        command = [
            self.codex_bin,
            "exec",
            "-",
            "--cd",
            str(working_dir),
            "--output-schema",
            str(schema_path),
            "-o",
            str(last_message_path),
        ]
        resolved_model = model if model is not None else self.model
        resolved_profile = profile if profile is not None else self.profile
        resolved_effort = reasoning_effort if reasoning_effort is not None else self.reasoning_effort
        resolved_config_overrides = list(self.config_overrides)
        if config_overrides:
            resolved_config_overrides.extend(config_overrides)
        if resolved_model:
            command.extend(["-m", resolved_model])
        if resolved_profile:
            command.extend(["--profile", resolved_profile])
        if resolved_effort:
            command.extend(["-c", f'model_reasoning_effort="{resolved_effort}"'])
        for override in resolved_config_overrides:
            command.extend(["-c", override])
        if self.full_auto:
            command.append("--full-auto")

        started = time.monotonic()
        heartbeat_stop = threading.Event()
        heartbeat_thread = None
        if progress_callback is not None and heartbeat_seconds > 0:
            def heartbeat_loop() -> None:
                while not heartbeat_stop.wait(heartbeat_seconds):
                    progress_callback(
                        "heartbeat",
                        {
                            "command": command,
                            "artifact_dir": str(artifact_dir),
                            "elapsed_seconds": round(time.monotonic() - started, 1),
                        },
                    )

            heartbeat_thread = threading.Thread(target=heartbeat_loop, daemon=True)
            heartbeat_thread.start()
        try:
            completed = None
            attempt_records: list[dict[str, object]] = []
            retries = max(0, int(max_retries))
            total_attempts = retries + 1
            for attempt_index in range(total_attempts):
                if last_message_path.exists():
                    last_message_path.unlink()
                attempt_started = time.monotonic()
                if progress_callback is not None and attempt_index > 0:
                    progress_callback(
                        "retrying",
                        {
                            "command": command,
                            "artifact_dir": str(artifact_dir),
                            "attempt": attempt_index + 1,
                            "max_attempts": total_attempts,
                        },
                    )
                completed = subprocess.run(
                    command,
                    cwd=working_dir,
                    input=prompt,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )
                raw_message = ""
                if last_message_path.exists():
                    raw_message = last_message_path.read_text(encoding="utf-8").strip()
                if not raw_message:
                    raw_message = completed.stdout.strip()
                retryable = self._is_retryable_disconnect(
                    returncode=completed.returncode,
                    stderr=completed.stderr,
                    raw_message=raw_message,
                )
                attempt_records.append(
                    {
                        "attempt": attempt_index + 1,
                        "returncode": completed.returncode,
                        "duration_seconds": round(time.monotonic() - attempt_started, 3),
                        "retryable_disconnect": retryable,
                        "stdout_tail": completed.stdout.splitlines()[-20:],
                        "stderr_tail": completed.stderr.splitlines()[-40:],
                    }
                )
                if not retryable or attempt_index == total_attempts - 1:
                    break
                time.sleep(min(2 ** attempt_index, 5))
        finally:
            heartbeat_stop.set()
            if heartbeat_thread is not None:
                heartbeat_thread.join(timeout=1.0)
        duration = round(time.monotonic() - started, 3)
        assert completed is not None
        write_json(attempts_path, {"attempts": attempt_records})

        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")

        raw_message = ""
        if last_message_path.exists():
            raw_message = last_message_path.read_text(encoding="utf-8").strip()
        if not raw_message:
            raw_message = completed.stdout.strip()

        parsed_message: dict
        try:
            parsed_message = json.loads(raw_message)
        except json.JSONDecodeError as exc:
            parsed_message = {
                "stage": "unparsed",
                "summary": "codex output did not parse as JSON",
                "selected_topics": [],
                "skills_touched": [],
                "resources_touched": [],
                "tests_run": [],
                "slurm_jobs": [],
                "repo_changes": [],
                "blockers": [str(exc)],
                "next_steps": [],
                "raw_message": raw_message,
            }

        return CodexInvocationResult(
            command=command,
            returncode=completed.returncode,
            duration_seconds=duration,
            stdout=completed.stdout,
            stderr=completed.stderr,
            raw_message=raw_message,
            parsed_message=parsed_message,
            prompt_path=prompt_path,
            schema_path=schema_path,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            last_message_path=last_message_path,
        )

    @staticmethod
    def _is_retryable_disconnect(*, returncode: int, stderr: str, raw_message: str) -> bool:
        stderr_lower = stderr.lower()
        raw_lower = raw_message.lower()
        transient_markers = [
            "stream disconnected before completion",
            "channel closed",
            "error sending request for url",
            "failed to queue rollout items",
        ]
        has_transient_marker = any(marker in stderr_lower for marker in transient_markers)
        no_useful_message = not raw_lower or "codex output did not parse as json" in raw_lower
        return has_transient_marker and no_useful_message
