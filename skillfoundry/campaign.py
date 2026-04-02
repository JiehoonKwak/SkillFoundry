"""Long-running checkpointable framework campaigns."""

from __future__ import annotations

import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path

from .orchestrator import SkillAutomationFramework
from .repository import write_json


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class FrameworkCampaignRunner:
    """Run checkpointable discovery + evaluation campaigns over the framework."""

    def __init__(self, framework: SkillAutomationFramework) -> None:
        self.framework = framework
        self.repo = framework.repo

    def _campaign_dir(self, label: str) -> Path:
        path = self.framework.state_dir / "campaigns" / label
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _state_paths(self, label: str) -> dict[str, Path]:
        base = self._campaign_dir(label)
        return {
            "base": base,
            "state": base / "campaign_state.json",
            "status": base / "status.json",
            "summary": base / "summary.md",
            "log": base / "campaign.log",
        }

    def _logger(self, label: str, log_path: Path) -> logging.Logger:
        logger = logging.getLogger(f"sciskill_framework.campaign.{label}")
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

    def _current_skill_slugs(self) -> list[str]:
        return sorted(skill["slug"] for skill in self.repo.load_skills())

    def _existing_campaign_run_dirs(self, label: str, kind: str) -> list[Path]:
        pattern = f"*-{label}-{kind}-*"
        return sorted(path for path in self.framework.reports_dir.glob(pattern) if path.is_dir())

    def _load_or_init_state(
        self,
        *,
        label: str,
        focus_terms: list[str],
        cycle_focus_limit: int,
        cycle_stage_workers: int,
        background_validation_limit: int,
        background_validation_workers: int,
        evaluation_batch_size: int,
        evaluation_workers: int,
        layer1_fix_attempts: int,
        layer2_optimize_attempts: int,
        verification_mode: str,
        full_verification_every: int,
        max_iterations: int,
        max_runtime_minutes: int,
        stop_buffer_minutes: int,
        extra_context: str | None,
    ) -> dict:
        paths = self._state_paths(label)
        if paths["state"].exists():
            state = json.loads(paths["state"].read_text(encoding="utf-8"))
            state["resumed_at"] = _utc_now()
            return state
        baseline_skills = self._current_skill_slugs()
        return {
            "label": label,
            "status": "initialized",
            "started_at": _utc_now(),
            "updated_at": _utc_now(),
            "focus_terms": list(focus_terms),
            "config": {
                "cycle_focus_limit": cycle_focus_limit,
                "cycle_stage_workers": cycle_stage_workers,
                "background_validation_limit": background_validation_limit,
                "background_validation_workers": background_validation_workers,
                "evaluation_batch_size": evaluation_batch_size,
                "evaluation_workers": evaluation_workers,
                "layer1_fix_attempts": layer1_fix_attempts,
                "layer2_optimize_attempts": layer2_optimize_attempts,
                "verification_mode": verification_mode,
                "full_verification_every": full_verification_every,
                "max_iterations": max_iterations,
                "max_runtime_minutes": max_runtime_minutes,
                "stop_buffer_minutes": stop_buffer_minutes,
                "extra_context": extra_context,
            },
            "baseline_skill_slugs": baseline_skills,
            "designed_skill_slugs": [],
            "new_skill_slugs": [],
            "completed_evaluations": {},
            "cycle_runs": [],
            "evaluation_runs": [],
            "periodic_verifications": [],
            "errors": [],
        }

    def _reconcile_state_from_artifacts(self, state: dict, logger: logging.Logger) -> None:
        tracked_cycle_dirs = {entry.get("run_dir") for entry in state.get("cycle_runs", [])}
        for run_dir in self._existing_campaign_run_dirs(state["label"], "cycle"):
            manifest_path = run_dir / "manifest.json"
            if not manifest_path.exists() or str(run_dir) in tracked_cycle_dirs:
                continue
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            self._record_cycle(state, payload)
            tracked_cycle_dirs.add(str(run_dir))
            logger.info("Reconciled existing cycle run %s into campaign state", run_dir.name)

        tracked_eval_dirs = {entry.get("run_dir") for entry in state.get("evaluation_runs", [])}
        for run_dir in self._existing_campaign_run_dirs(state["label"], "eval"):
            manifest_path = run_dir / "manifest.json"
            if not manifest_path.exists() or str(run_dir) in tracked_eval_dirs:
                continue
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            self._record_evaluation(state, payload)
            tracked_eval_dirs.add(str(run_dir))
            logger.info("Reconciled existing evaluation run %s into campaign state", run_dir.name)

        incomplete_resume = None
        completed = set(state.get("completed_evaluations", {}).keys())
        for run_dir in self._existing_campaign_run_dirs(state["label"], "eval"):
            if (run_dir / "manifest.json").exists():
                continue
            skill_root = run_dir / "skills"
            if not skill_root.exists():
                continue
            skill_slugs = sorted(path.name for path in skill_root.iterdir() if path.is_dir())
            if not skill_slugs:
                continue
            if all(slug in completed for slug in skill_slugs):
                continue
            incomplete_resume = {
                "run_label": run_dir.name,
                "run_dir": str(run_dir),
                "skill_slugs": skill_slugs,
            }
            logger.info(
                "Discovered incomplete evaluation run %s with %s skills; it will be resumed before the next cycle",
                run_dir.name,
                len(skill_slugs),
            )
            break
        if incomplete_resume is not None:
            state["resume_evaluation"] = incomplete_resume
        else:
            state.pop("resume_evaluation", None)

    def _pending_evaluation_slugs(self, state: dict) -> list[str]:
        completed = set(state.get("completed_evaluations", {}).keys())
        return [slug for slug in self._current_skill_slugs() if slug not in completed]

    def _latest_summary(self, state: dict) -> dict:
        current_skills = self._current_skill_slugs()
        completed = state.get("completed_evaluations", {})
        pending = self._pending_evaluation_slugs(state)
        summary = self.framework.status(
            focus_limit=state["config"]["cycle_focus_limit"],
            focus_terms=state.get("focus_terms"),
        )
        active_step_name = None
        active_step_summary = None
        active_run_label = None
        active_run_path = self.framework.state_dir / "latest_active_run.json"
        if active_run_path.exists():
            try:
                active_run = json.loads(active_run_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                active_run = {}
            if active_run.get("active"):
                current_step = active_run.get("current_step") or {}
                active_step_name = current_step.get("name")
                active_step_summary = current_step.get("summary")
                active_run_label = active_run.get("run_label")
        return {
            "label": state["label"],
            "status": state["status"],
            "started_at": state.get("started_at"),
            "updated_at": _utc_now(),
            "focus_terms": state.get("focus_terms", []),
            "iterations_completed": len(state.get("cycle_runs", [])),
            "cycle_runs_completed": len(state.get("cycle_runs", [])),
            "evaluation_runs_completed": len(state.get("evaluation_runs", [])),
            "baseline_skill_count": len(state.get("baseline_skill_slugs", [])),
            "current_skill_count": len(current_skills),
            "designed_skill_count": len(state.get("designed_skill_slugs", [])),
            "new_skill_count": len(state.get("new_skill_slugs", [])),
            "evaluated_skill_count": len(completed),
            "layer1_pass_count": sum(1 for result in completed.values() if result.get("layer1_passed")),
            "layer2_pass_count": sum(1 for result in completed.values() if result.get("layer2_passed")),
            "pending_evaluation_count": len(pending),
            "pending_evaluation_preview": pending[:10],
            "last_cycle_run_dir": state.get("cycle_runs", [{}])[-1].get("run_dir") if state.get("cycle_runs") else None,
            "last_evaluation_run_dir": state.get("evaluation_runs", [{}])[-1].get("run_dir") if state.get("evaluation_runs") else None,
            "active_run_label": active_run_label,
            "active_step_name": active_step_name,
            "active_step_summary": active_step_summary,
            "tree_summary": summary["summary"],
            "focus_leaves": summary["focus_leaves"],
        }

    def _write_state(self, state: dict) -> dict:
        paths = self._state_paths(state["label"])
        state["updated_at"] = _utc_now()
        write_json(paths["state"], state)
        status_payload = self._latest_summary(state)
        write_json(paths["status"], status_payload)
        self._write_summary_markdown(paths["summary"], state, status_payload)
        return status_payload

    def _write_summary_markdown(self, path: Path, state: dict, status_payload: dict) -> None:
        lines = [
            f"# Campaign Summary: {state['label']}",
            "",
            "## Aggregate",
            "",
            "| Metric | Value |",
            "| --- | --- |",
            f"| Status | {status_payload['status']} |",
            f"| Focus terms | {', '.join(status_payload['focus_terms']) or 'none'} |",
            f"| Cycle runs completed | {status_payload['cycle_runs_completed']} |",
            f"| Evaluation runs completed | {status_payload['evaluation_runs_completed']} |",
            f"| Current skills in registry | {status_payload['current_skill_count']} |",
            f"| Designed or updated skills | {status_payload['designed_skill_count']} |",
            f"| New skills | {status_payload['new_skill_count']} |",
            f"| Evaluated skills | {status_payload['evaluated_skill_count']} |",
            f"| Layer 1 passed | {status_payload['layer1_pass_count']} |",
            f"| Layer 2 passed | {status_payload['layer2_pass_count']} |",
            f"| Pending evaluation skills | {status_payload['pending_evaluation_count']} |",
            f"| Active framework run | {status_payload.get('active_run_label') or 'n/a'} |",
            f"| Active step | {status_payload.get('active_step_name') or 'n/a'} |",
            "",
            "## Latest Tree Snapshot",
            "",
            "| Metric | Value |",
            "| --- | --- |",
            f"| Resources | {status_payload['tree_summary']['resource_count']} |",
            f"| Skills | {status_payload['tree_summary']['skill_count']} |",
            f"| Covered leaves | {status_payload['tree_summary']['covered_leaf_count']} |",
            f"| Frontier leaves | {status_payload['tree_summary']['frontier_leaf_count']} |",
            f"| TODO leaves | {status_payload['tree_summary']['todo_leaf_count']} |",
            "",
            "## Campaign Skill Table",
            "",
            "| Skill | New | Designed/Updated | Layer 1 | Layer 2 | Last evaluation run |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        completed = state.get("completed_evaluations", {})
        new_skills = set(state.get("new_skill_slugs", []))
        designed = set(state.get("designed_skill_slugs", []))
        tracked = sorted(set(completed) | new_skills | designed)
        for slug in tracked:
            result = completed.get(slug, {})
            lines.append(
                f"| {slug} | {'yes' if slug in new_skills else 'no'} | "
                f"{'yes' if slug in designed else 'no'} | "
                f"{'pass' if result.get('layer1_passed') else ('fail' if slug in completed else 'n/a')} | "
                f"{'pass' if result.get('layer2_passed') else ('fail' if slug in completed else 'n/a')} | "
                f"{result.get('evaluation_run_label', 'n/a')} |"
            )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _record_cycle(self, state: dict, payload: dict) -> None:
        designed = {
            skill
            for record in payload.get("stage_results", [])
            for skill in record.get("codex", {}).get("parsed_message", {}).get("skills_touched", [])
        }
        state["designed_skill_slugs"] = sorted(set(state.get("designed_skill_slugs", [])) | designed)
        current_skills = set(self._current_skill_slugs())
        baseline = set(state.get("baseline_skill_slugs", []))
        state["new_skill_slugs"] = sorted(current_skills - baseline)
        cycle_record = {
            "run_label": payload["run_label"],
            "run_dir": payload["run_dir"],
            "verification_failures": sum(1 for record in payload.get("verification", []) if record.get("returncode") != 0),
            "background_validation_failed": (payload.get("background_validation") or {}).get("failed_skill_count", 0),
            "focus_terms": payload.get("focus_terms", []),
            "final_summary": payload.get("final_status", {}).get("summary", {}),
        }
        cycle_runs = state.setdefault("cycle_runs", [])
        for index, existing in enumerate(cycle_runs):
            if existing.get("run_dir") == payload["run_dir"]:
                cycle_runs[index] = cycle_record
                break
        else:
            cycle_runs.append(cycle_record)

    def _record_evaluation(self, state: dict, payload: dict) -> None:
        completed = state.setdefault("completed_evaluations", {})
        layer1_by_slug = {record["skill_slug"]: record for record in payload.get("layer1", [])}
        layer2_by_slug = {record["skill_slug"]: record for record in payload.get("layer2", [])}
        novelty_by_slug = {record["skill_slug"]: record for record in payload.get("novelty", [])}
        for skill in payload.get("selected_skills", []):
            slug = skill["slug"]
            completed[slug] = {
                "skill_id": skill["skill_id"],
                "domain": skill["domain"],
                "status": skill["status"],
                "layer1_passed": bool(layer1_by_slug.get(slug, {}).get("passed")),
                "layer2_passed": bool(layer2_by_slug.get(slug, {}).get("meets_expectation")),
                "novelty_checked": slug in novelty_by_slug,
                "evaluation_run_label": payload["run_label"],
                "evaluation_run_dir": payload["run_dir"],
            }
        eval_record = {
            "run_label": payload["run_label"],
            "run_dir": payload["run_dir"],
            "selected_skill_count": payload.get("evaluation_summary", {}).get("selected_skill_count", 0),
            "layer1_passed": payload.get("evaluation_summary", {}).get("layer1_passed", 0),
            "layer2_advantage_passed": payload.get("evaluation_summary", {}).get("layer2_advantage_passed", 0),
            "novelty_checked": payload.get("evaluation_summary", {}).get("novelty_checked", 0),
        }
        evaluation_runs = state.setdefault("evaluation_runs", [])
        for index, existing in enumerate(evaluation_runs):
            if existing.get("run_dir") == payload["run_dir"]:
                evaluation_runs[index] = eval_record
                break
        else:
            evaluation_runs.append(eval_record)

    def _run_periodic_verification(self, state: dict, iteration_index: int) -> None:
        full_every = max(0, int(state["config"].get("full_verification_every", 0)))
        if full_every <= 0 or iteration_index % full_every != 0:
            return
        run_dir = self._campaign_dir(state["label"]) / "periodic-verification" / f"iteration-{iteration_index:03d}"
        records = self.repo.execute_verification("full", run_dir=run_dir, timeout=self.framework.verification_timeout)
        state.setdefault("periodic_verifications", []).append(
            {
                "iteration": iteration_index,
                "run_dir": str(run_dir),
                "records": records,
            }
        )

    def run_campaign(
        self,
        *,
        label: str,
        focus_terms: list[str],
        cycle_focus_limit: int = 12,
        cycle_stage_workers: int = 6,
        background_validation_limit: int = 24,
        background_validation_workers: int = 6,
        evaluation_batch_size: int = 24,
        evaluation_workers: int = 6,
        layer1_fix_attempts: int = 1,
        layer2_optimize_attempts: int = 1,
        verification_mode: str = "standard",
        full_verification_every: int = 0,
        max_iterations: int = 1000,
        max_runtime_minutes: int = 450,
        stop_buffer_minutes: int = 20,
        extra_context: str | None = None,
    ) -> dict:
        paths = self._state_paths(label)
        logger = self._logger(label, paths["log"])
        state = self._load_or_init_state(
            label=label,
            focus_terms=focus_terms,
            cycle_focus_limit=cycle_focus_limit,
            cycle_stage_workers=cycle_stage_workers,
            background_validation_limit=background_validation_limit,
            background_validation_workers=background_validation_workers,
            evaluation_batch_size=evaluation_batch_size,
            evaluation_workers=evaluation_workers,
            layer1_fix_attempts=layer1_fix_attempts,
            layer2_optimize_attempts=layer2_optimize_attempts,
            verification_mode=verification_mode,
            full_verification_every=full_verification_every,
            max_iterations=max_iterations,
            max_runtime_minutes=max_runtime_minutes,
            stop_buffer_minutes=stop_buffer_minutes,
            extra_context=extra_context,
        )
        self._reconcile_state_from_artifacts(state, logger)
        state["status"] = "running"
        self._write_state(state)

        started = time.monotonic()
        deadline = started + max(1, max_runtime_minutes) * 60
        stop_buffer_seconds = max(0, stop_buffer_minutes) * 60

        logger.info("Campaign %s started with focus terms: %s", label, ", ".join(focus_terms) or "none")
        logger.info(
            "Workers: cycle=%s background_validation=%s evaluation=%s",
            cycle_stage_workers,
            background_validation_workers,
            evaluation_workers,
        )

        try:
            while state.get("resume_evaluation") or len(state.get("cycle_runs", [])) < max_iterations:
                remaining = deadline - time.monotonic()
                if remaining <= stop_buffer_seconds:
                    logger.info("Stopping before walltime buffer: %.1f minutes remaining", remaining / 60.0)
                    state["status"] = "checkpointed"
                    break

                resume_evaluation = state.get("resume_evaluation")
                if resume_evaluation:
                    run_dir = Path(resume_evaluation["run_dir"])
                    logger.info(
                        "Resuming incomplete evaluation run %s for %s skills",
                        resume_evaluation["run_label"],
                        len(resume_evaluation["skill_slugs"]),
                    )
                    evaluation_payload = self.framework.evaluate_skills(
                        skill_slugs=resume_evaluation["skill_slugs"],
                        include_all=False,
                        limit=len(resume_evaluation["skill_slugs"]),
                        workers=evaluation_workers,
                        verification_mode="none",
                        run_label=resume_evaluation["run_label"],
                        run_dir=run_dir,
                        layer1_fix_attempts=layer1_fix_attempts,
                        layer2_optimize_attempts=layer2_optimize_attempts,
                        novelty_check=True,
                        extra_context=extra_context,
                    )
                    self._record_evaluation(state, evaluation_payload)
                    state.pop("resume_evaluation", None)
                    status_payload = self._write_state(state)
                    logger.info(
                        "Resumed evaluation %s finished: layer1=%s layer2=%s pending=%s",
                        resume_evaluation["run_label"],
                        evaluation_payload["evaluation_summary"]["layer1_passed"],
                        evaluation_payload["evaluation_summary"]["layer2_advantage_passed"],
                        status_payload["pending_evaluation_count"],
                    )
                    continue

                iteration_index = len(state.get("cycle_runs", [])) + 1
                logger.info("Iteration %s starting", iteration_index)
                before_new_count = len(state.get("new_skill_slugs", []))
                before_designed_count = len(state.get("designed_skill_slugs", []))
                cycle_payload = self.framework.run_cycle(
                    loops=1,
                    stages=["tree_check", "resource_search", "skill_build", "skill_test", "refresh"],
                    focus_limit=cycle_focus_limit,
                    focus_terms=focus_terms,
                    stage_workers=cycle_stage_workers,
                    background_validation_limit=background_validation_limit,
                    background_validation_workers=background_validation_workers,
                    verification_mode=verification_mode,
                    run_label=f"{label}-cycle-{iteration_index:03d}",
                    extra_context=extra_context,
                )
                self._record_cycle(state, cycle_payload)
                self._write_state(state)
                cycle_added_activity = (
                    len(state.get("new_skill_slugs", [])) > before_new_count
                    or len(state.get("designed_skill_slugs", [])) > before_designed_count
                )
                logger.info(
                    "Iteration %s cycle finished: run_dir=%s new_skills=%s designed=%s",
                    iteration_index,
                    cycle_payload["run_dir"],
                    len(state.get("new_skill_slugs", [])),
                    len(state.get("designed_skill_slugs", [])),
                )

                pending = self._pending_evaluation_slugs(state)
                if pending:
                    batch = pending[: max(1, evaluation_batch_size)]
                    logger.info("Evaluating %s pending skills in batch %s", len(batch), iteration_index)
                    evaluation_payload = self.framework.evaluate_skills(
                        skill_slugs=batch,
                        include_all=False,
                        limit=len(batch),
                        workers=evaluation_workers,
                        verification_mode="none",
                        run_label=f"{label}-eval-{iteration_index:03d}",
                        layer1_fix_attempts=layer1_fix_attempts,
                        layer2_optimize_attempts=layer2_optimize_attempts,
                        novelty_check=True,
                        extra_context=extra_context,
                    )
                    self._record_evaluation(state, evaluation_payload)
                    logger.info(
                        "Evaluation batch %s finished: layer1_passed=%s layer2_passed=%s",
                        iteration_index,
                        evaluation_payload["evaluation_summary"]["layer1_passed"],
                        evaluation_payload["evaluation_summary"]["layer2_advantage_passed"],
                    )
                else:
                    logger.info("No pending evaluation skills remain")

                self._run_periodic_verification(state, iteration_index)
                status_payload = self._write_state(state)
                logger.info(
                    "Checkpoint written: designed=%s new=%s layer1=%s layer2=%s pending=%s",
                    status_payload["designed_skill_count"],
                    status_payload["new_skill_count"],
                    status_payload["layer1_pass_count"],
                    status_payload["layer2_pass_count"],
                    status_payload["pending_evaluation_count"],
                )

                if status_payload["pending_evaluation_count"] == 0 and not cycle_added_activity:
                    state["status"] = "completed"
                    break
                if status_payload["pending_evaluation_count"] == 0 and remaining <= stop_buffer_seconds + 300:
                    state["status"] = "completed"
                    break

            if state["status"] == "running":
                state["status"] = "completed" if self._pending_evaluation_slugs(state) == [] else "checkpointed"
        except Exception as exc:
            logger.exception("Campaign failed: %s", exc)
            state.setdefault("errors", []).append({"at": _utc_now(), "error": str(exc)})
            state["status"] = "failed"
            self._write_state(state)
            raise

        status_payload = self._write_state(state)
        logger.info("Campaign %s finished with status=%s", label, state["status"])
        return {
            "mode": "campaign",
            "label": label,
            "campaign_dir": str(paths["base"]),
            "state_path": str(paths["state"]),
            "status_path": str(paths["status"]),
            "summary_path": str(paths["summary"]),
            "log_path": str(paths["log"]),
            "campaign_summary": status_payload,
        }

    def campaign_status(self, *, label: str) -> dict:
        paths = self._state_paths(label)
        if not paths["state"].exists():
            raise FileNotFoundError(f"No campaign state found for label: {label}")
        state = json.loads(paths["state"].read_text(encoding="utf-8"))
        logger = self._logger(label, paths["log"])
        self._reconcile_state_from_artifacts(state, logger)
        status_payload = self._write_state(state)
        return {
            "mode": "campaign_status",
            "label": label,
            "campaign_dir": str(paths["base"]),
            "state_path": str(paths["state"]),
            "status_path": str(paths["status"]),
            "summary_path": str(paths["summary"]),
            "log_path": str(paths["log"]),
            "campaign_summary": status_payload,
        }
