"""Main orchestration entrypoint for the codex-exec automation framework."""

from __future__ import annotations

import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import re

from .codex_exec import CodexExecRunner
from .evaluation import HierarchicalSkillEvaluator
from .live_state import FrameworkRunStateWriter
from .prompts import (
    DEFAULT_CYCLE_STAGES,
    PARALLEL_LEAF_STAGES,
    build_design_skill_prompt,
    build_parallel_leaf_stage_prompt,
    build_stage_prompt,
    stage_output_schema,
)
from .repository import SciSkillRepository, write_json

DEFAULT_STAGE_MODEL_ROUTING: dict[str, dict[str, str]] = {
    "tree_check": {"model": "gpt-5.4", "reasoning_effort": "medium"},
    "resource_search": {"model": "gpt-5.4", "reasoning_effort": "high"},
    "skill_build": {"model": "gpt-5.4", "reasoning_effort": "high"},
    "skill_test": {"model": "gpt-5.4", "reasoning_effort": "high"},
    "refresh": {"model": "gpt-5.4", "reasoning_effort": "medium"},
    "design_skill": {"model": "gpt-5.4", "reasoning_effort": "high"},
    "layer1_fix": {"model": "gpt-5.4", "reasoning_effort": "medium"},
    "layer2_benchmark": {"model": "gpt-5.4-mini", "reasoning_effort": "medium"},
    "layer2_optimize": {"model": "gpt-5.4", "reasoning_effort": "medium"},
    "novelty_check": {"model": "gpt-5.4-mini", "reasoning_effort": "medium"},
}


class SkillAutomationFramework:
    """Automate the resource-search -> skill-build -> test loop via `codex exec`."""

    def __init__(
        self,
        repo_root: str | Path,
        *,
        codex_bin: str = "codex",
        state_dir: str | Path | None = None,
        reports_dir: str | Path | None = None,
        workspace_root: str | Path | None = None,
        model: str | None = None,
        reasoning_effort: str | None = None,
        profile: str | None = None,
        full_auto: bool = False,
        codex_config_overrides: list[str] | None = None,
        stage_timeout: int = 3600,
        verification_timeout: int = 3600,
    ) -> None:
        self.repo = SciSkillRepository(Path(repo_root))
        self.state_dir = Path(state_dir) if state_dir is not None else self.repo.scratch / "framework"
        self.reports_dir = Path(reports_dir) if reports_dir is not None else self.repo.reports / "framework-runs"
        self.workspace_root = Path(workspace_root) if workspace_root is not None else self.state_dir / "workspaces"
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        self.workspace_root.mkdir(parents=True, exist_ok=True)
        self.stage_timeout = stage_timeout
        self.verification_timeout = verification_timeout
        self.default_codex_options = {
            "model": model,
            "reasoning_effort": reasoning_effort,
            "profile": profile,
        }
        self.runner = CodexExecRunner(
            codex_bin=codex_bin,
            model=model,
            reasoning_effort=reasoning_effort,
            profile=profile,
            full_auto=full_auto,
            config_overrides=codex_config_overrides,
        )
        self.live_state = FrameworkRunStateWriter(self.state_dir)
        self.evaluator = HierarchicalSkillEvaluator(
            repo=self.repo,
            runner=self.runner,
            state_dir=self.state_dir,
            live_state=self.live_state,
            stage_timeout=self.stage_timeout,
            verification_timeout=self.verification_timeout,
            codex_options_for_stage=self._codex_options_for_stage,
        )

    def _timestamp(self) -> str:
        return datetime.now().strftime("%Y%m%d-%H%M%S")

    def _run_dir(self, prefix: str) -> Path:
        run_dir = self.reports_dir / f"{self._timestamp()}-{prefix}"
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir

    def _write_state(self, name: str, payload: dict) -> Path:
        path = self.state_dir / name
        write_json(path, payload)
        return path

    def _verification_progress_callback(self, event: str, record: dict) -> None:
        step_name = f"verify:{record['label']}"
        if event == "start":
            self.live_state.mark_step_running(
                step_name,
                kind="verification",
                summary=f"Running {record['label']} verification",
                detail={"command": " ".join(record["command"])},
            )
            self.live_state.append_log(
                kind="command",
                label=step_name,
                text=" ".join(record["command"]),
            )
            return
        self.live_state.append_lines(
            kind="stdout",
            label=step_name,
            lines=record.get("stdout_tail", []),
        )
        self.live_state.append_lines(
            kind="stderr",
            label=step_name,
            lines=record.get("stderr_tail", []),
        )
        self.live_state.mark_step_complete(
            step_name,
            summary=f"{record['label']} verification finished with code {record['returncode']}",
            health="clean" if record.get("returncode") == 0 else "attention",
            blockers=[] if record.get("returncode") == 0 else [f"{record['label']} verification failed"],
        )

    def _codex_progress_callback(self, label: str):
        def callback(event: str, record: dict) -> None:
            if event != "heartbeat":
                return
            elapsed = int(record.get("elapsed_seconds", 0))
            self.live_state.append_log(
                kind="info",
                label=label,
                text=f"Still running after {elapsed}s",
            )

        return callback

    def _slugify(self, value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-") or "worker"

    def _codex_options_for_stage(self, stage: str) -> dict[str, str]:
        options = dict(DEFAULT_STAGE_MODEL_ROUTING.get(stage, {}))
        for key, value in self.default_codex_options.items():
            if value is not None:
                options[key] = value
        return {key: value for key, value in options.items() if value is not None}

    def _test_name_safe(self, name: str) -> bool:
        return bool(re.fullmatch(r"test_[a-z0-9_]+\.py", name))

    def _worker_workspace_path(self, run_dir: Path, focus_leaf) -> Path:
        return self.workspace_root / run_dir.name / f"{focus_leaf.domain_slug}--{self._slugify(focus_leaf.topic_slug)}"

    def _workspace_ignore(self, dirpath: str, names: list[str]) -> set[str]:
        ignored = {name for name in names if name == "__pycache__" or name.endswith(".pyc")}
        try:
            rel_dir = Path(dirpath).resolve().relative_to(self.repo.root)
        except ValueError:
            rel_dir = Path(".")
        if rel_dir == Path("."):
            ignored.update({name for name in names if name in {".git", "scratch", "reports", "output", ".nextflow", ".nextflow.log", ".nextflow.log.1"}})
        if rel_dir == Path("slurm"):
            ignored.update({name for name in names if name in {"envs", "logs"}})
        return ignored

    def _ensure_worker_workspace(self, run_dir: Path, focus_leaf) -> Path:
        workspace = self._worker_workspace_path(run_dir, focus_leaf)
        if workspace.exists():
            return workspace
        workspace.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            self.repo.root,
            workspace,
            ignore=self._workspace_ignore,
        )
        env_target = self.repo.root / "slurm" / "envs"
        env_link = workspace / "slurm" / "envs"
        env_link.parent.mkdir(parents=True, exist_ok=True)
        if env_target.exists() and not env_link.exists():
            env_link.symlink_to(env_target, target_is_directory=True)
        return workspace

    def _is_worker_sync_path(self, rel_path: Path) -> bool:
        rel = rel_path.as_posix()
        if rel.startswith("skills/"):
            return True
        if rel.startswith(("tests/smoke/", "tests/regression/", "tests/integration/", "tests/slurm/")) and self._test_name_safe(rel_path.name):
            return True
        if rel.startswith("slurm/jobs/"):
            return True
        return False

    def _workspace_changed_paths(self, workspace: Path) -> list[Path]:
        changed: list[Path] = []
        for prefix in ("skills", "tests", "slurm/jobs"):
            base = workspace / prefix
            if not base.exists():
                continue
            for path in base.rglob("*"):
                if not path.is_file() or path.is_symlink():
                    continue
                rel_path = path.relative_to(workspace)
                if not self._is_worker_sync_path(rel_path):
                    continue
                main_path = self.repo.root / rel_path
                if not main_path.exists() or path.read_bytes() != main_path.read_bytes():
                    changed.append(rel_path)
        deduped = sorted({path.as_posix(): path for path in changed}.values(), key=lambda item: item.as_posix())
        return deduped

    def _sync_worker_workspaces(self, run_dir: Path, focus_leaves: list, *, workspace_records: dict[str, Path]) -> dict:
        sync_dir = run_dir / "workspace-sync"
        sync_dir.mkdir(parents=True, exist_ok=True)
        ownership: dict[str, list[str]] = {}
        worker_changes: dict[str, list[str]] = {}
        for focus_leaf in focus_leaves:
            workspace = workspace_records.get(focus_leaf.topic_slug)
            if workspace is None or not workspace.exists():
                continue
            changes = [path.as_posix() for path in self._workspace_changed_paths(workspace)]
            worker_changes[focus_leaf.topic_slug] = changes
            for rel_path in changes:
                ownership.setdefault(rel_path, []).append(focus_leaf.topic_slug)

        conflict_paths = sorted(path for path, owners in ownership.items() if len(owners) > 1)
        synced_paths: list[str] = []
        for rel_path, owners in sorted(ownership.items()):
            if len(owners) > 1:
                self.live_state.append_log(
                    kind="warning",
                    label="workspace_sync",
                    text=f"Skipped conflicted worker path {rel_path} touched by {', '.join(sorted(owners))}",
                )
                continue
            owner = owners[0]
            workspace = workspace_records[owner]
            src = workspace / rel_path
            dest = self.repo.root / rel_path
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            synced_paths.append(rel_path)

        record = {
            "workspace_root": str(self.workspace_root / run_dir.name),
            "worker_changes": worker_changes,
            "synced_paths": synced_paths,
            "conflict_paths": conflict_paths,
        }
        write_json(sync_dir / "sync_record.json", record)
        self.live_state.append_log(
            kind="success" if not conflict_paths else "warning",
            label="workspace_sync",
            text=f"Synced {len(synced_paths)} worker-owned path(s) with {len(conflict_paths)} conflict(s).",
        )
        return record

    def _background_validation_progress_callback(self, lane_name: str, skill_slug: str):
        label = f"{lane_name}:{skill_slug}"

        def callback(event: str, record: dict) -> None:
            if event == "start":
                self.live_state.append_log(
                    kind="command",
                    label=label,
                    text=" ".join(record["command"]),
                )
                return
            self.live_state.append_lines(kind="stdout", label=label, lines=record.get("stdout_tail", []))
            self.live_state.append_lines(kind="stderr", label=label, lines=record.get("stderr_tail", []))
            self.live_state.append_log(
                kind="success" if record.get("returncode") == 0 else "error",
                label=label,
                text=(
                    f"{record.get('target_label') or 'skill check'} finished with code "
                    f"{record.get('returncode')} in {record.get('duration_seconds', 0.0):.3f}s"
                ),
            )

        return callback

    def _run_background_validation_lane(
        self,
        *,
        run_dir: Path,
        skill_limit: int,
        workers: int,
    ) -> dict:
        selected_skills = self.repo.select_skills(limit=skill_limit)
        smoke_map = self.repo.smoke_map(selected_skills)
        validation_dir = run_dir / "background-validation"
        validation_dir.mkdir(parents=True, exist_ok=True)
        self.live_state.mark_step_running(
            "parallel_validation",
            kind="stage",
            summary=f"Running background validation for {len(selected_skills)} skills",
        )

        def validate_one(position_skill):
            position, skill = position_skill
            self.live_state.append_log(
                kind="info",
                label="parallel_validation",
                text=f"[{position}/{len(selected_skills)}] Validating {skill['slug']}",
            )
            record = self.repo.execute_skill_check(
                skill,
                artifact_dir=validation_dir / skill["slug"],
                timeout=self.verification_timeout,
                smoke_map=smoke_map,
                progress_callback=self._background_validation_progress_callback("parallel_validation", skill["slug"]),
            )
            return position, record

        results: list[tuple[int, dict]] = []
        with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
            futures = [executor.submit(validate_one, item) for item in enumerate(selected_skills, start=1)]
            for future in as_completed(futures):
                results.append(future.result())

        ordered = [record for _, record in sorted(results, key=lambda item: item[0])]
        failures = [record["slug"] for record in ordered if record.get("returncode") != 0]
        summary = {
            "selected_skill_count": len(selected_skills),
            "workers": max(1, workers),
            "failed_skill_count": len(failures),
            "failed_skills": failures,
            "records": ordered,
        }
        write_json(validation_dir / "summary.json", summary)
        self.live_state.mark_step_complete(
            "parallel_validation",
            summary=f"Background validation finished for {len(selected_skills)} skills; {len(failures)} failed.",
            health="clean" if not failures else "attention",
            blockers=[f"{slug} failed background validation" for slug in failures],
        )
        return summary

    def _execute_parallel_leaf_stage(
        self,
        *,
        stage: str,
        stage_index: int,
        focus_leaves: list,
        run_dir: Path,
        extra_context: str | None,
        workers: int,
        workspace_records: dict[str, Path],
    ) -> dict:
        pre_summary = self.repo.summary(focus_limit=len(focus_leaves), focus_terms=[leaf.topic_slug for leaf in focus_leaves])
        artifact_dir = run_dir / "stages" / f"{stage_index:02d}-{stage}"
        workers_dir = artifact_dir / "workers"
        workers_dir.mkdir(parents=True, exist_ok=True)

        def run_one(position_leaf):
            position, focus_leaf = position_leaf
            worker_label = f"{stage}:{focus_leaf.topic_slug}"
            workspace = self._ensure_worker_workspace(run_dir, focus_leaf)
            workspace_records[focus_leaf.topic_slug] = workspace
            self.live_state.append_log(
                kind="info",
                label=stage,
                text=f"[{position}/{len(focus_leaves)}] Worker assigned to {focus_leaf.topic_slug}",
            )
            worker_artifact_dir = workers_dir / f"{position:02d}-{self._slugify(focus_leaf.topic_slug)}"
            prompt = build_parallel_leaf_stage_prompt(
                stage=stage,
                repo_summary=pre_summary,
                focus_leaf=focus_leaf,
                artifact_dir=str(worker_artifact_dir),
                extra_context=extra_context,
            )
            result = self.runner.run(
                prompt=prompt,
                output_schema=stage_output_schema(stage),
                working_dir=workspace,
                artifact_dir=worker_artifact_dir,
                **self._codex_options_for_stage(stage),
                timeout=self.stage_timeout,
                progress_callback=self._codex_progress_callback(worker_label),
            )
            result_dict = result.to_dict()
            self.live_state.append_lines(kind="stdout", label=worker_label, lines=result_dict.get("stdout_tail", []))
            self.live_state.append_lines(kind="stderr", label=worker_label, lines=result_dict.get("stderr_tail", []))
            self.live_state.append_log(
                kind="success" if result_dict.get("returncode") == 0 else "error",
                label=worker_label,
                text=result_dict["parsed_message"].get("summary", f"{stage} worker completed"),
            )
            return position, {
                "worker_index": position,
                "focus_leaf": focus_leaf.to_dict(),
                "workspace": str(workspace),
                "candidate_sync_paths": [path.as_posix() for path in self._workspace_changed_paths(workspace)],
                "codex": result_dict,
            }

        collected: list[tuple[int, dict]] = []
        with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
            futures = [executor.submit(run_one, item) for item in enumerate(focus_leaves, start=1)]
            for future in as_completed(futures):
                collected.append(future.result())

        worker_results = [record for _, record in sorted(collected, key=lambda item: item[0])]
        aggregate = {
            "stage": stage,
            "summary": (
                f"Parallel {stage} completed for {len(worker_results)} focus leaves with "
                f"{sum(1 for record in worker_results if record['codex'].get('returncode') != 0 or record['codex']['parsed_message'].get('blockers'))} worker(s) needing follow-up."
            ),
            "selected_topics": [record["focus_leaf"]["topic_slug"] for record in worker_results],
            "skills_touched": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("skills_touched", [])}),
            "resources_touched": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("resources_touched", [])}),
            "tests_run": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("tests_run", [])}),
            "slurm_jobs": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("slurm_jobs", [])}),
            "repo_changes": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("repo_changes", [])}),
            "blockers": [
                f"{record['focus_leaf']['topic_slug']}: {blocker}"
                for record in worker_results
                for blocker in record["codex"]["parsed_message"].get("blockers", [])
            ],
            "next_steps": sorted({item for record in worker_results for item in record["codex"]["parsed_message"].get("next_steps", [])}),
        }
        stage_record = {
            "stage": stage,
            "index": stage_index,
            "parallel": True,
            "worker_count_requested": max(1, workers),
            "pre_summary": pre_summary.to_dict(),
            "post_summary": pre_summary.to_dict(),
            "focus_leaves": [leaf.to_dict() for leaf in focus_leaves],
            "worker_results": worker_results,
            "codex": {
                "command": [],
                "returncode": 0 if not aggregate["blockers"] else 1,
                "duration_seconds": 0.0,
                "stdout_tail": [],
                "stderr_tail": [],
                "prompt_path": None,
                "schema_path": None,
                "last_message_path": None,
                "parsed_message": aggregate,
            },
        }
        write_json(artifact_dir / "result.json", stage_record)
        return stage_record

    def status(self, focus_limit: int = 8, focus_terms: list[str] | None = None) -> dict:
        summary = self.repo.summary(focus_limit=focus_limit, focus_terms=focus_terms)
        payload = {
            "summary": summary.to_dict(),
            "focus_leaves": [leaf.to_dict() for leaf in summary.focus_leaves],
        }
        self._write_state("latest_status.json", payload)
        return payload

    def _execute_stage(
        self,
        *,
        stage: str,
        stage_index: int,
        focus_limit: int,
        focus_terms: list[str] | None,
        run_dir: Path,
        extra_context: str | None = None,
    ) -> dict:
        pre_summary = self.repo.summary(focus_limit=focus_limit, focus_terms=focus_terms)
        artifact_dir = run_dir / "stages" / f"{stage_index:02d}-{stage}"
        prompt = build_stage_prompt(
            stage=stage,
            repo_summary=pre_summary,
            focus_leaves=pre_summary.focus_leaves,
            extra_context=extra_context,
        )
        result = self.runner.run(
            prompt=prompt,
            output_schema=stage_output_schema(stage),
            working_dir=self.repo.root,
            artifact_dir=artifact_dir,
            **self._codex_options_for_stage(stage),
            timeout=self.stage_timeout,
            progress_callback=self._codex_progress_callback(stage),
        )
        post_summary = self.repo.summary(focus_limit=focus_limit, focus_terms=focus_terms)
        stage_record = {
            "stage": stage,
            "index": stage_index,
            "pre_summary": pre_summary.to_dict(),
            "post_summary": post_summary.to_dict(),
            "focus_leaves": [leaf.to_dict() for leaf in pre_summary.focus_leaves],
            "codex": result.to_dict(),
        }
        write_json(artifact_dir / "result.json", stage_record)
        return stage_record

    def run_cycle(
        self,
        *,
        loops: int = 1,
        stages: list[str] | tuple[str, ...] | None = None,
        focus_limit: int = 8,
        focus_terms: list[str] | None = None,
        stage_workers: int = 1,
        background_validation_limit: int = 0,
        background_validation_workers: int = 1,
        verification_mode: str = "standard",
        run_label: str = "cycle",
        extra_context: str | None = None,
    ) -> dict:
        stage_list = list(stages or DEFAULT_CYCLE_STAGES)
        run_dir = self._run_dir(run_label)
        verification_steps = [{"name": f"verify:{label}", "kind": "verification", "status": "pending", "summary": ""} for label, _ in self.repo.verification_commands(verification_mode)]
        extra_steps = []
        if background_validation_limit > 0:
            extra_steps.append({"name": "parallel_validation", "kind": "stage", "status": "pending", "summary": ""})
        self.live_state.start_run(
            mode="cycle",
            run_label=run_label,
            run_dir=str(run_dir),
            verification_mode=verification_mode,
            steps=[{"name": stage, "kind": "stage", "status": "pending", "summary": ""} for stage in stage_list] + extra_steps + verification_steps,
        )
        manifest = {
            "mode": "cycle",
            "run_label": run_label,
            "run_dir": str(run_dir),
            "loops": loops,
            "stages_requested": stage_list,
            "focus_terms": list(focus_terms or []),
            "stage_workers": max(1, stage_workers),
            "background_validation_limit": background_validation_limit,
            "background_validation_workers": max(1, background_validation_workers),
            "verification_mode": verification_mode,
            "initial_status": self.status(focus_limit=focus_limit, focus_terms=focus_terms),
            "stage_results": [],
            "background_validation": None,
            "workspace_sync": None,
            "verification": [],
        }
        write_json(run_dir / "initial_status.json", manifest["initial_status"])

        background_executor = ThreadPoolExecutor(max_workers=1) if background_validation_limit > 0 else None
        background_future = None
        try:
            stage_index = 1
            workspace_records: dict[str, Path] = {}
            sync_focus_leaves: list = []
            pending_workspace_sync = False
            for loop_index in range(loops):
                for stage in stage_list:
                    if stage == "refresh" and pending_workspace_sync and workspace_records:
                        manifest["workspace_sync"] = self._sync_worker_workspaces(
                            run_dir,
                            sync_focus_leaves,
                            workspace_records=workspace_records,
                        )
                        pending_workspace_sync = False
                    self.live_state.mark_step_running(
                        stage,
                        kind="stage",
                        summary=f"Loop {loop_index + 1} / {loops} running {stage}",
                        detail={"loop_index": loop_index + 1},
                    )
                    self.live_state.append_log(
                        kind="info",
                        label=stage,
                        text=f"Starting {stage} (loop {loop_index + 1} of {loops})",
                    )
                    if background_executor is not None and background_future is None and stage != "tree_check":
                        background_future = background_executor.submit(
                            self._run_background_validation_lane,
                            run_dir=run_dir,
                            skill_limit=background_validation_limit,
                            workers=background_validation_workers,
                        )
                    if stage in PARALLEL_LEAF_STAGES and max(1, stage_workers) > 1:
                        focus_leaves = self.repo.summary(focus_limit=focus_limit, focus_terms=focus_terms).focus_leaves
                        sync_focus_leaves = focus_leaves
                        stage_record = self._execute_parallel_leaf_stage(
                            stage=stage,
                            stage_index=stage_index,
                            focus_leaves=focus_leaves,
                            run_dir=run_dir,
                            extra_context=f"Loop {loop_index + 1} of {loops}. {extra_context or ''}".strip(),
                            workers=stage_workers,
                            workspace_records=workspace_records,
                        )
                        pending_workspace_sync = True
                    else:
                        stage_record = self._execute_stage(
                            stage=stage,
                            stage_index=stage_index,
                            focus_limit=focus_limit,
                            focus_terms=focus_terms,
                            run_dir=run_dir,
                            extra_context=f"Loop {loop_index + 1} of {loops}. {extra_context or ''}".strip(),
                        )
                    stage_record["loop_index"] = loop_index + 1
                    manifest["stage_results"].append(stage_record)
                    parsed = stage_record["codex"]["parsed_message"]
                    self.live_state.append_lines(
                        kind="stdout",
                        label=stage,
                        lines=stage_record["codex"].get("stdout_tail", []),
                    )
                    self.live_state.append_lines(
                        kind="stderr",
                        label=stage,
                        lines=stage_record["codex"].get("stderr_tail", []),
                    )
                    self.live_state.mark_step_complete(
                        stage,
                        summary=parsed.get("summary", f"{stage} completed"),
                        health="clean" if not parsed.get("blockers") and stage_record["codex"].get("returncode") == 0 else "attention",
                        blockers=parsed.get("blockers", []),
                    )
                    stage_index += 1

            if pending_workspace_sync and workspace_records:
                manifest["workspace_sync"] = self._sync_worker_workspaces(
                    run_dir,
                    sync_focus_leaves,
                    workspace_records=workspace_records,
                )

            if background_future is not None:
                manifest["background_validation"] = background_future.result()
            if background_executor is not None:
                background_executor.shutdown(wait=True)
                background_executor = None

            manifest["verification"] = self.repo.execute_verification(
                verification_mode,
                run_dir=run_dir,
                timeout=self.verification_timeout,
                progress_callback=self._verification_progress_callback,
            )
            manifest["final_status"] = self.status(focus_limit=focus_limit, focus_terms=focus_terms)
            write_json(run_dir / "manifest.json", manifest)
            self._write_state("latest_run.json", manifest)
            verification_failures = sum(1 for record in manifest["verification"] if record.get("returncode") != 0)
            self.live_state.finish(
                status="completed" if verification_failures == 0 else "failed",
                summary=f"Cycle finished with {len(manifest['stage_results'])} stages and {verification_failures} verification failure(s).",
            )
            return manifest
        except Exception as exc:
            self.live_state.finish(status="failed", summary=f"Cycle aborted: {exc}")
            raise
        finally:
            if background_future is not None and not background_future.done():
                background_future.cancel()
            if background_executor is not None:
                background_executor.shutdown(wait=False, cancel_futures=True)

    def design_skill(
        self,
        *,
        task_prompt: str,
        focus_limit: int = 8,
        focus_terms: list[str] | None = None,
        verification_mode: str = "standard",
        run_label: str = "design-skill",
        extra_context: str | None = None,
    ) -> dict:
        run_dir = self._run_dir(run_label)
        pre_summary = self.repo.summary(focus_limit=focus_limit, focus_terms=focus_terms)
        verification_steps = [{"name": f"verify:{label}", "kind": "verification", "status": "pending", "summary": ""} for label, _ in self.repo.verification_commands(verification_mode)]
        self.live_state.start_run(
            mode="design_skill",
            run_label=run_label,
            run_dir=str(run_dir),
            verification_mode=verification_mode,
            steps=[{"name": "design_skill", "kind": "stage", "status": "pending", "summary": ""}] + verification_steps,
        )
        prompt = build_design_skill_prompt(
            task_prompt=task_prompt,
            repo_summary=pre_summary,
            focus_leaves=pre_summary.focus_leaves,
            extra_context=extra_context,
        )
        try:
            self.live_state.mark_step_running(
                "design_skill",
                kind="stage",
                summary="Designing a task-specific skill path",
                detail={"task_prompt": task_prompt},
            )
            result = self.runner.run(
                prompt=prompt,
                output_schema=stage_output_schema("design_skill"),
                working_dir=self.repo.root,
                artifact_dir=run_dir / "stages" / "01-design-skill",
                **self._codex_options_for_stage("design_skill"),
                timeout=self.stage_timeout,
                progress_callback=self._codex_progress_callback("design_skill"),
            )
            result_dict = result.to_dict()
            self.live_state.append_lines(kind="stdout", label="design_skill", lines=result_dict.get("stdout_tail", []))
            self.live_state.append_lines(kind="stderr", label="design_skill", lines=result_dict.get("stderr_tail", []))
            self.live_state.mark_step_complete(
                "design_skill",
                summary=result_dict["parsed_message"].get("summary", "design_skill completed"),
                health="clean" if not result_dict["parsed_message"].get("blockers") and result_dict.get("returncode") == 0 else "attention",
                blockers=result_dict["parsed_message"].get("blockers", []),
            )
            post_summary = self.repo.summary(focus_limit=focus_limit)
            verification = self.repo.execute_verification(
                verification_mode,
                run_dir=run_dir,
                timeout=self.verification_timeout,
                progress_callback=self._verification_progress_callback,
            )
            manifest = {
                "mode": "design_skill",
                "run_label": run_label,
                "run_dir": str(run_dir),
                "task_prompt": task_prompt,
                "focus_terms": list(focus_terms or []),
                "verification_mode": verification_mode,
                "initial_status": pre_summary.to_dict(),
                "stage_result": {
                    "stage": "design_skill",
                    "pre_summary": pre_summary.to_dict(),
                    "post_summary": post_summary.to_dict(),
                    "focus_leaves": [leaf.to_dict() for leaf in pre_summary.focus_leaves],
                    "codex": result_dict,
                },
                "verification": verification,
                "final_status": self.status(focus_limit=focus_limit, focus_terms=focus_terms),
            }
            write_json(run_dir / "manifest.json", manifest)
            self._write_state("latest_design_skill.json", manifest)
            verification_failures = sum(1 for record in verification if record.get("returncode") != 0)
            self.live_state.finish(
                status="completed" if verification_failures == 0 else "failed",
                summary=f"Design-skill finished with {verification_failures} verification failure(s).",
            )
            return manifest
        except Exception as exc:
            self.live_state.finish(status="failed", summary=f"Design-skill aborted: {exc}")
            raise

    def evaluate_skills(
        self,
        *,
        skill_slugs: list[str] | None = None,
        include_all: bool = False,
        limit: int = 8,
        workers: int = 1,
        verification_mode: str = "standard",
        run_label: str = "evaluate-skills",
        run_dir: Path | None = None,
        layer1_fix_attempts: int = 1,
        layer2_optimize_attempts: int = 1,
        novelty_check: bool = True,
        extra_context: str | None = None,
    ) -> dict:
        selected_run_dir = run_dir if run_dir is not None else self._run_dir(run_label)
        payload = self.evaluator.evaluate_skills(
            run_dir=selected_run_dir,
            skill_slugs=skill_slugs,
            include_all=include_all,
            limit=limit,
            workers=workers,
            verification_mode=verification_mode,
            layer1_fix_attempts=layer1_fix_attempts,
            layer2_optimize_attempts=layer2_optimize_attempts,
            novelty_check=novelty_check,
            extra_context=extra_context,
        )
        self._write_state("latest_evaluation.json", payload)
        return payload
