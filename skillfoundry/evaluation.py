"""Hierarchical evaluation flow for existing or newly developed skills."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .codex_exec import CodexExecRunner
from .live_state import FrameworkRunStateWriter
from .prompts import (
    build_layer1_fix_prompt,
    build_layer2_benchmark_prompt,
    build_layer2_optimize_prompt,
    build_novelty_check_prompt,
    layer1_fix_output_schema,
    layer2_benchmark_output_schema,
    layer2_optimize_output_schema,
    novelty_check_output_schema,
)
from .repository import SciSkillRepository, write_json


class HierarchicalSkillEvaluator:
    """Run layered correctness, benchmark, and novelty evaluation for skills."""

    def __init__(
        self,
        *,
        repo: SciSkillRepository,
        runner: CodexExecRunner,
        state_dir: Path,
        live_state: FrameworkRunStateWriter | None = None,
        stage_timeout: int = 3600,
        skill_check_timeout: int = 1800,
        verification_timeout: int = 3600,
        codex_options_for_stage=None,
    ) -> None:
        self.repo = repo
        self.runner = runner
        self.state_dir = Path(state_dir)
        self.stage_timeout = stage_timeout
        self.skill_check_timeout = skill_check_timeout
        self.verification_timeout = verification_timeout
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.live_state = live_state or FrameworkRunStateWriter(self.state_dir)
        self.codex_options_for_stage = codex_options_for_stage or (lambda stage: {})

    def _stage_summary(
        self,
        *,
        stage: str,
        summary: str,
        selected_topics: list[str],
        skills_touched: list[str],
        tests_run: list[str],
        blockers: list[str],
        next_steps: list[str],
        resources_touched: list[str] | None = None,
        repo_changes: list[str] | None = None,
    ) -> dict:
        return {
            "stage": stage,
            "codex": {
                "returncode": 0 if not blockers else 1,
                "parsed_message": {
                    "stage": stage,
                    "summary": summary,
                    "selected_topics": selected_topics,
                    "skills_touched": skills_touched,
                    "resources_touched": resources_touched or [],
                    "tests_run": tests_run,
                    "slurm_jobs": [],
                    "repo_changes": repo_changes or [],
                    "blockers": blockers,
                    "next_steps": next_steps,
                },
            },
        }

    def _run_codex(
        self,
        *,
        stage_name: str,
        prompt: str,
        output_schema: dict,
        artifact_dir: Path,
        label: str,
    ) -> dict:
        def progress_callback(event: str, record: dict) -> None:
            if event != "heartbeat":
                return
            elapsed = int(record.get("elapsed_seconds", 0))
            self.live_state.append_log(
                kind="info",
                label=label,
                text=f"Still running after {elapsed}s",
            )

        result = self.runner.run(
            prompt=prompt,
            output_schema=output_schema,
            working_dir=self.repo.root,
            artifact_dir=artifact_dir,
            **self.codex_options_for_stage(stage_name),
            timeout=self.stage_timeout,
            progress_callback=progress_callback,
        )
        return result.to_dict()

    def _layer_record_path(self, *, run_dir: Path, skill_slug: str, layer: str) -> Path:
        filenames = {
            "layer1": "layer1_record.json",
            "layer2": "layer2_record.json",
            "novelty": "novelty_record.json",
        }
        return run_dir / "skills" / skill_slug / layer / filenames[layer]

    def _load_existing_record(self, *, run_dir: Path, skill_slug: str, layer: str) -> dict | None:
        path = self._layer_record_path(run_dir=run_dir, skill_slug=skill_slug, layer=layer)
        if not path.exists():
            return None
        try:
            return self.repo.load_json(path)
        except Exception:
            return None

    def _layer1_passed(self, record: dict) -> bool:
        final_check = record.get("final_check", {})
        return final_check.get("returncode") == 0

    def _benchmark_passed(self, record: dict) -> bool:
        final_benchmark = record.get("final_benchmark", {})
        if "parsed_message" in final_benchmark:
            return final_benchmark.get("returncode") == 0 and bool(final_benchmark.get("parsed_message", {}).get("meets_expectation"))
        codex = final_benchmark.get("codex", {})
        parsed = codex.get("parsed_message", {})
        return codex.get("returncode") == 0 and bool(parsed.get("meets_expectation"))

    def _skill_check_progress_callback(self, *, stage_name: str, skill_slug: str):
        label = f"{stage_name}:{skill_slug}"

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
            level = "success" if record.get("returncode") == 0 else "error"
            self.live_state.append_log(
                kind=level,
                label=label,
                text=(
                    f"{record.get('target_label') or 'skill check'} finished with code "
                    f"{record.get('returncode')} in {record.get('duration_seconds', 0.0):.3f}s"
                ),
            )

        return callback

    def _run_layer1(
        self,
        *,
        skill: dict,
        run_dir: Path,
        repo_summary,
        smoke_map: dict[str, list[str]],
        fix_attempts: int,
        extra_context: str | None,
        position: int,
        total: int,
    ) -> dict:
        skill_dir = run_dir / "skills" / skill["slug"] / "layer1"
        self.live_state.append_log(
            kind="info",
            label="layer1_correctness",
            text=f"[{position}/{total}] Running correctness check for {skill['slug']}",
        )
        initial_check = self.repo.execute_skill_check(
            skill,
            artifact_dir=skill_dir / "check-00",
            timeout=self.skill_check_timeout,
            smoke_map=smoke_map,
            progress_callback=self._skill_check_progress_callback(
                stage_name="layer1",
                skill_slug=skill["slug"],
            ),
        )
        record = {
            "skill_id": skill["skill_id"],
            "skill_slug": skill["slug"],
            "skill_path": skill["path"],
            "initial_check": initial_check,
            "fix_attempts": [],
            "final_check": initial_check,
            "passed": initial_check.get("returncode") == 0,
        }
        current_check = initial_check
        for attempt in range(1, fix_attempts + 1):
            if current_check.get("returncode") == 0:
                break
            fix_dir = skill_dir / f"fix-{attempt:02d}"
            fix_prompt = build_layer1_fix_prompt(
                skill=skill,
                failure_record=current_check,
                repo_summary=repo_summary,
                artifact_dir=str(fix_dir),
                extra_context=extra_context,
            )
            fix_result = self._run_codex(
                stage_name="layer1_fix",
                prompt=fix_prompt,
                output_schema=layer1_fix_output_schema(),
                artifact_dir=fix_dir,
                label=f"layer1_fix:{skill['slug']}",
            )
            self.live_state.append_lines(kind="stdout", label=skill["slug"], lines=fix_result.get("stdout_tail", []))
            self.live_state.append_lines(kind="stderr", label=skill["slug"], lines=fix_result.get("stderr_tail", []))
            retry_check = self.repo.execute_skill_check(
                skill,
                artifact_dir=skill_dir / f"check-{attempt:02d}",
                timeout=self.skill_check_timeout,
                smoke_map=smoke_map,
                progress_callback=self._skill_check_progress_callback(
                    stage_name=f"layer1-recheck-{attempt:02d}",
                    skill_slug=skill["slug"],
                ),
            )
            record["fix_attempts"].append(
                {
                    "attempt": attempt,
                    "fix": fix_result,
                    "recheck": retry_check,
                }
            )
            current_check = retry_check

        record["final_check"] = current_check
        record["passed"] = current_check.get("returncode") == 0
        write_json(skill_dir / "layer1_record.json", record)
        return record

    def _run_layer2(
        self,
        *,
        skill: dict,
        run_dir: Path,
        repo_summary,
        smoke_map: dict[str, list[str]],
        optimize_attempts: int,
        extra_context: str | None,
        position: int,
        total: int,
    ) -> dict:
        skill_dir = run_dir / "skills" / skill["slug"] / "layer2"
        self.live_state.append_log(
            kind="info",
            label="layer2_advantage",
            text=f"[{position}/{total}] Benchmarking {skill['slug']} against a no-skill baseline",
        )
        benchmark_dir = skill_dir / "benchmark-00"
        benchmark_prompt = build_layer2_benchmark_prompt(
            skill=skill,
            repo_summary=repo_summary,
            artifact_dir=str(benchmark_dir),
            extra_context=extra_context,
        )
        benchmark_result = self._run_codex(
            stage_name="layer2_benchmark",
            prompt=benchmark_prompt,
            output_schema=layer2_benchmark_output_schema(),
            artifact_dir=benchmark_dir,
            label=f"layer2_benchmark:{skill['slug']}",
        )
        self.live_state.append_lines(kind="stdout", label=skill["slug"], lines=benchmark_result.get("stdout_tail", []))
        self.live_state.append_lines(kind="stderr", label=skill["slug"], lines=benchmark_result.get("stderr_tail", []))
        record = {
            "skill_id": skill["skill_id"],
            "skill_slug": skill["slug"],
            "skill_path": skill["path"],
            "initial_benchmark": benchmark_result,
            "optimize_attempts": [],
            "final_benchmark": benchmark_result,
            "meets_expectation": self._benchmark_passed({"final_benchmark": benchmark_result}),
        }
        current_benchmark = benchmark_result

        for attempt in range(1, optimize_attempts + 1):
            parsed = current_benchmark.get("parsed_message", {})
            if current_benchmark.get("returncode") == 0 and parsed.get("meets_expectation"):
                break
            optimize_dir = skill_dir / f"optimize-{attempt:02d}"
            optimize_prompt = build_layer2_optimize_prompt(
                skill=skill,
                benchmark_result=current_benchmark,
                repo_summary=repo_summary,
                artifact_dir=str(optimize_dir),
                extra_context=extra_context,
            )
            optimize_result = self._run_codex(
                stage_name="layer2_optimize",
                prompt=optimize_prompt,
                output_schema=layer2_optimize_output_schema(),
                artifact_dir=optimize_dir,
                label=f"layer2_optimize:{skill['slug']}",
            )
            self.live_state.append_lines(kind="stdout", label=skill["slug"], lines=optimize_result.get("stdout_tail", []))
            self.live_state.append_lines(kind="stderr", label=skill["slug"], lines=optimize_result.get("stderr_tail", []))
            post_opt_check = self.repo.execute_skill_check(
                skill,
                artifact_dir=skill_dir / f"post-opt-check-{attempt:02d}",
                timeout=self.skill_check_timeout,
                smoke_map=smoke_map,
                progress_callback=self._skill_check_progress_callback(
                    stage_name=f"layer2-post-opt-{attempt:02d}",
                    skill_slug=skill["slug"],
                ),
            )
            rerun_benchmark = None
            if post_opt_check.get("returncode") == 0:
                rerun_dir = skill_dir / f"benchmark-{attempt:02d}"
                rerun_prompt = build_layer2_benchmark_prompt(
                    skill=skill,
                    repo_summary=repo_summary,
                    artifact_dir=str(rerun_dir),
                    extra_context=extra_context,
                )
                rerun_benchmark = self._run_codex(
                    stage_name="layer2_benchmark",
                    prompt=rerun_prompt,
                    output_schema=layer2_benchmark_output_schema(),
                    artifact_dir=rerun_dir,
                    label=f"layer2_benchmark:{skill['slug']}",
                )
                self.live_state.append_lines(kind="stdout", label=skill["slug"], lines=rerun_benchmark.get("stdout_tail", []))
                self.live_state.append_lines(kind="stderr", label=skill["slug"], lines=rerun_benchmark.get("stderr_tail", []))
                current_benchmark = rerun_benchmark
            record["optimize_attempts"].append(
                {
                    "attempt": attempt,
                    "optimize": optimize_result,
                    "post_opt_check": post_opt_check,
                    "benchmark": rerun_benchmark,
                }
            )
            if post_opt_check.get("returncode") != 0:
                break

        record["final_benchmark"] = current_benchmark
        record["meets_expectation"] = self._benchmark_passed(record)
        write_json(skill_dir / "layer2_record.json", record)
        return record

    def _run_novelty(
        self,
        *,
        skill: dict,
        run_dir: Path,
        repo_summary,
        extra_context: str | None,
        position: int,
        total: int,
    ) -> dict:
        skill_dir = run_dir / "skills" / skill["slug"] / "novelty"
        candidates = self.repo.local_similarity_candidates(skill, limit=5)
        self.live_state.append_log(
            kind="info",
            label="novelty_check",
            text=f"[{position}/{total}] Checking novelty and overlap for {skill['slug']}",
        )
        prompt = build_novelty_check_prompt(
            skill=skill,
            local_candidates=candidates,
            repo_summary=repo_summary,
            artifact_dir=str(skill_dir),
            extra_context=extra_context,
        )
        novelty_result = self._run_codex(
            stage_name="novelty_check",
            prompt=prompt,
            output_schema=novelty_check_output_schema(),
            artifact_dir=skill_dir,
            label=f"novelty_check:{skill['slug']}",
        )
        self.live_state.append_lines(kind="stdout", label=skill["slug"], lines=novelty_result.get("stdout_tail", []))
        self.live_state.append_lines(kind="stderr", label=skill["slug"], lines=novelty_result.get("stderr_tail", []))
        record = {
            "skill_id": skill["skill_id"],
            "skill_slug": skill["slug"],
            "skill_path": skill["path"],
            "local_similarity_candidates": candidates,
            "novelty": novelty_result,
        }
        write_json(skill_dir / "novelty_record.json", record)
        return record

    def evaluate_skills(
        self,
        *,
        run_dir: Path,
        skill_slugs: list[str] | None = None,
        include_all: bool = False,
        limit: int = 8,
        workers: int = 1,
        verification_mode: str = "standard",
        layer1_fix_attempts: int = 1,
        layer2_optimize_attempts: int = 1,
        novelty_check: bool = True,
        extra_context: str | None = None,
    ) -> dict:
        selected_skills = self.repo.select_skills(
            skill_slugs=skill_slugs,
            include_all=include_all,
            limit=limit,
        )
        run_dir.mkdir(parents=True, exist_ok=True)
        repo_summary = self.repo.summary(focus_limit=min(max(limit, 1), 16))
        smoke_map = self.repo.smoke_map(selected_skills)

        manifest = {
            "mode": "evaluate_skills",
            "run_label": run_dir.name,
            "run_dir": str(run_dir),
            "selected_skills": [
                {
                    "skill_id": skill["skill_id"],
                    "slug": skill["slug"],
                    "domain": skill["domain"],
                    "status": skill["status"],
                    "path": skill["path"],
                }
                for skill in selected_skills
            ],
            "evaluation_config": {
                "include_all": include_all,
                "limit": limit,
                "workers": max(1, workers),
                "verification_mode": verification_mode,
                "layer1_fix_attempts": layer1_fix_attempts,
                "layer2_optimize_attempts": layer2_optimize_attempts,
                "novelty_check": novelty_check,
            },
            "initial_status": {
                "summary": repo_summary.to_dict(),
                "focus_leaves": [leaf.to_dict() for leaf in repo_summary.focus_leaves],
            },
            "layer1": [],
            "layer2": [],
            "novelty": [],
            "human_expert_evaluation": {
                "implemented": False,
                "note": "Human expert evaluation is intentionally left as a manual layer outside this automation framework.",
            },
            "verification": [],
            "stage_results": [],
        }
        verification_steps = [{"name": f"verify:{label}", "kind": "verification", "status": "pending", "summary": ""} for label, _ in self.repo.verification_commands(verification_mode)]
        self.live_state.start_run(
            mode="evaluate_skills",
            run_label=run_dir.name,
            run_dir=str(run_dir),
            verification_mode=verification_mode,
            steps=[
                {"name": "layer1_correctness", "kind": "stage", "status": "pending", "summary": ""},
                {"name": "layer2_advantage", "kind": "stage", "status": "pending", "summary": ""},
                {"name": "novelty_check", "kind": "stage", "status": "pending", "summary": ""},
            ] + verification_steps,
            selected_skills=[skill["slug"] for skill in selected_skills],
        )
        if len(selected_skills) > 20:
            self.live_state.append_log(
                kind="warning",
                label="evaluate_skills",
                text=(
                    f"Large evaluation requested for {len(selected_skills)} skills. "
                    "This is a long-running job because layer1, layer2, and novelty "
                        "review can each touch every selected skill."
                ),
            )
        worker_count = max(1, workers)

        try:
            self.live_state.mark_step_running(
                "layer1_correctness",
                kind="stage",
                summary="Running sandbox / Slurm-facing correctness checks",
            )
            existing_layer1 = {
                skill["slug"]: record
                for skill in selected_skills
                if (record := self._load_existing_record(run_dir=run_dir, skill_slug=skill["slug"], layer="layer1")) is not None
            }
            missing_layer1_skills = [skill for skill in selected_skills if skill["slug"] not in existing_layer1]
            if existing_layer1:
                self.live_state.append_log(
                    kind="info",
                    label="layer1_correctness",
                    text=f"Reused {len(existing_layer1)} existing layer1 record(s) from {run_dir.name}.",
                )
            layer1_results = list(existing_layer1.values())
            if missing_layer1_skills:
                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    futures = [
                        executor.submit(
                            self._run_layer1,
                            skill=skill,
                            run_dir=run_dir,
                            repo_summary=repo_summary,
                            smoke_map=smoke_map,
                            fix_attempts=layer1_fix_attempts,
                            extra_context=extra_context,
                            position=position,
                            total=len(selected_skills),
                        )
                        for position, skill in enumerate(missing_layer1_skills, start=1)
                    ]
                    layer1_results.extend(future.result() for future in as_completed(futures))
            layer1_by_slug = {record["skill_slug"]: record for record in layer1_results}
            manifest["layer1"] = [layer1_by_slug[skill["slug"]] for skill in selected_skills if skill["slug"] in layer1_by_slug]
            layer1_failed = [record["skill_slug"] for record in manifest["layer1"] if not record["passed"]]
            self.live_state.mark_step_complete(
                "layer1_correctness",
                summary=f"Layer 1 completed for {len(manifest['layer1'])} skills; {len(layer1_failed)} still failing.",
                health="clean" if not layer1_failed else "attention",
                blockers=[f"{slug} failed correctness" for slug in layer1_failed],
            )

            layer1_passed_skills = [
                skill
                for skill in selected_skills
                if next(
                    (
                        record["passed"]
                        for record in manifest["layer1"]
                        if record["skill_slug"] == skill["slug"]
                    ),
                    False,
                )
            ]
            self.live_state.mark_step_running(
                "layer2_advantage",
                kind="stage",
                summary="Running skill-vs-baseline benchmark comparisons",
            )
            existing_layer2 = {
                skill["slug"]: record
                for skill in layer1_passed_skills
                if (record := self._load_existing_record(run_dir=run_dir, skill_slug=skill["slug"], layer="layer2")) is not None
            }
            missing_layer2_skills = [skill for skill in layer1_passed_skills if skill["slug"] not in existing_layer2]
            if existing_layer2:
                self.live_state.append_log(
                    kind="info",
                    label="layer2_advantage",
                    text=f"Reused {len(existing_layer2)} existing layer2 record(s) from {run_dir.name}.",
                )
            layer2_results = list(existing_layer2.values())
            if missing_layer2_skills:
                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    futures = [
                        executor.submit(
                            self._run_layer2,
                            skill=skill,
                            run_dir=run_dir,
                            repo_summary=repo_summary,
                            smoke_map=smoke_map,
                            optimize_attempts=layer2_optimize_attempts,
                            extra_context=extra_context,
                            position=position,
                            total=len(layer1_passed_skills),
                        )
                        for position, skill in enumerate(missing_layer2_skills, start=1)
                    ]
                    layer2_results.extend(future.result() for future in as_completed(futures))
            layer2_by_slug = {record["skill_slug"]: record for record in layer2_results}
            manifest["layer2"] = [layer2_by_slug[skill["slug"]] for skill in layer1_passed_skills if skill["slug"] in layer2_by_slug]
            layer2_failed = [record["skill_slug"] for record in manifest["layer2"] if not record["meets_expectation"]]
            self.live_state.mark_step_complete(
                "layer2_advantage",
                summary=f"Layer 2 benchmarked {len(manifest['layer2'])} skills; {len(layer2_failed)} still need optimization.",
                health="clean" if not layer2_failed else "attention",
                blockers=[f"{slug} did not beat the baseline" for slug in layer2_failed],
            )

            if novelty_check:
                self.live_state.mark_step_running(
                    "novelty_check",
                    kind="stage",
                    summary="Running overlap and novelty review",
                )
                existing_novelty = {
                    skill["slug"]: record
                    for skill in selected_skills
                    if (record := self._load_existing_record(run_dir=run_dir, skill_slug=skill["slug"], layer="novelty")) is not None
                }
                missing_novelty_skills = [skill for skill in selected_skills if skill["slug"] not in existing_novelty]
                if existing_novelty:
                    self.live_state.append_log(
                        kind="info",
                        label="novelty_check",
                        text=f"Reused {len(existing_novelty)} existing novelty record(s) from {run_dir.name}.",
                    )
                novelty_results = list(existing_novelty.values())
                if missing_novelty_skills:
                    with ThreadPoolExecutor(max_workers=worker_count) as executor:
                        futures = [
                            executor.submit(
                                self._run_novelty,
                                skill=skill,
                                run_dir=run_dir,
                                repo_summary=repo_summary,
                                extra_context=extra_context,
                                position=position,
                                total=len(selected_skills),
                            )
                            for position, skill in enumerate(missing_novelty_skills, start=1)
                        ]
                        novelty_results.extend(future.result() for future in as_completed(futures))
                novelty_by_slug = {record["skill_slug"]: record for record in novelty_results}
                manifest["novelty"] = [novelty_by_slug[skill["slug"]] for skill in selected_skills if skill["slug"] in novelty_by_slug]
                novelty_merge_review = [
                    record["skill_slug"]
                    for record in manifest["novelty"]
                    if record["novelty"]["parsed_message"].get("needs_merge_review")
                ]
                self.live_state.mark_step_complete(
                    "novelty_check",
                    summary=f"Novelty review completed for {len(manifest['novelty'])} skills; {len(novelty_merge_review)} need merge review.",
                    health="clean" if not novelty_merge_review else "attention",
                    blockers=[f"{slug} needs merge review" for slug in novelty_merge_review],
                )
            else:
                novelty_merge_review = []
                self.live_state.mark_step_complete(
                    "novelty_check",
                    summary="Novelty review skipped for this run.",
                    health="clean",
                    blockers=[],
                )

            layer1_failed = [record["skill_slug"] for record in manifest["layer1"] if not record["passed"]]
            layer2_failed = [record["skill_slug"] for record in manifest["layer2"] if not record["meets_expectation"]]
            manifest["evaluation_summary"] = {
                "selected_skill_count": len(selected_skills),
                "layer1_passed": len(manifest["layer1"]) - len(layer1_failed),
                "layer1_failed": len(layer1_failed),
                "layer1_failed_skills": layer1_failed,
                "layer2_evaluated": len(manifest["layer2"]),
                "layer2_advantage_passed": len(manifest["layer2"]) - len(layer2_failed),
                "layer2_advantage_failed": len(layer2_failed),
                "layer2_failed_skills": layer2_failed,
                "novelty_checked": len(manifest["novelty"]),
                "novelty_merge_review_count": len(novelty_merge_review),
                "novelty_merge_review_skills": novelty_merge_review,
            }

            layer1_tests = [
                record["final_check"]["target_label"]
                for record in manifest["layer1"]
                if record["final_check"].get("target_label")
            ]
            layer1_blockers = [f"{slug}: correctness check failed" for slug in layer1_failed]
            layer2_tests = []
            for record in manifest["layer2"]:
                parsed = record["final_benchmark"]["parsed_message"]
                layer2_tests.extend(parsed.get("tests_run", []))
            layer2_blockers = [f"{slug}: benchmark advantage not proven" for slug in layer2_failed]
            novelty_blockers = [f"{slug}: merge review recommended" for slug in novelty_merge_review]
            novelty_resources = []
            for record in manifest["novelty"]:
                novelty_resources.extend(record["novelty"]["parsed_message"].get("web_sources", []))

            manifest["stage_results"] = [
                self._stage_summary(
                    stage="layer1_correctness",
                    summary=(
                        f"Layer 1 checked {len(manifest['layer1'])} skills; "
                        f"{manifest['evaluation_summary']['layer1_passed']} passed correctness."
                    ),
                    selected_topics=[record["skill_slug"] for record in manifest["layer1"]],
                    skills_touched=[record["skill_slug"] for record in manifest["layer1"]],
                    tests_run=layer1_tests,
                    blockers=layer1_blockers,
                    next_steps=["Investigate failed correctness checks first."] if layer1_blockers else ["Proceed to benchmark validation."],
                ),
                self._stage_summary(
                    stage="layer2_advantage",
                    summary=(
                        f"Layer 2 benchmarked {len(manifest['layer2'])} skills; "
                        f"{manifest['evaluation_summary']['layer2_advantage_passed']} showed a measured advantage."
                    ),
                    selected_topics=[record["skill_slug"] for record in manifest["layer2"]],
                    skills_touched=[record["skill_slug"] for record in manifest["layer2"]],
                    tests_run=layer2_tests,
                    blockers=layer2_blockers,
                    next_steps=["Optimize skills that did not beat the baseline."] if layer2_blockers else ["Retain benchmark artifacts for reporting."],
                ),
                self._stage_summary(
                    stage="novelty_check",
                    summary=(
                        f"Novelty checks completed for {len(manifest['novelty'])} skills; "
                        f"{manifest['evaluation_summary']['novelty_merge_review_count']} need merge review."
                    ),
                    selected_topics=[record["skill_slug"] for record in manifest["novelty"]],
                    skills_touched=[record["skill_slug"] for record in manifest["novelty"]],
                    tests_run=[],
                    blockers=novelty_blockers,
                    next_steps=["Inspect overlap-heavy skills before creating more adjacent skills."] if novelty_blockers else ["Novelty surface looks acceptable for this batch."],
                    resources_touched=sorted(set(novelty_resources)),
                ),
            ]

            def verification_callback(event: str, record: dict) -> None:
                step_name = f"verify:{record['label']}"
                if event == "start":
                    self.live_state.mark_step_running(
                        step_name,
                        kind="verification",
                        summary=f"Running {record['label']} verification",
                        detail={"command": " ".join(record["command"])},
                    )
                    self.live_state.append_log(kind="command", label=step_name, text=" ".join(record["command"]))
                    return
                self.live_state.append_lines(kind="stdout", label=step_name, lines=record.get("stdout_tail", []))
                self.live_state.append_lines(kind="stderr", label=step_name, lines=record.get("stderr_tail", []))
                self.live_state.mark_step_complete(
                    step_name,
                    summary=f"{record['label']} verification finished with code {record['returncode']}",
                    health="clean" if record.get("returncode") == 0 else "attention",
                    blockers=[] if record.get("returncode") == 0 else [f"{record['label']} verification failed"],
                )

            manifest["verification"] = self.repo.execute_verification(
                verification_mode,
                run_dir=run_dir,
                timeout=self.verification_timeout,
                progress_callback=verification_callback,
            )
            manifest["final_status"] = {
                "summary": self.repo.summary(focus_limit=min(max(limit, 1), 16)).to_dict(),
            }

            write_json(run_dir / "manifest.json", manifest)
            write_json(self.state_dir / "latest_evaluation.json", manifest)
            verification_failures = sum(1 for record in manifest["verification"] if record.get("returncode") != 0)
            self.live_state.finish(
                status="completed" if verification_failures == 0 else "failed",
                summary=(
                    f"Hierarchical evaluation finished for {len(selected_skills)} skills with "
                    f"{verification_failures} verification failure(s)."
                ),
            )
            return manifest
        except Exception as exc:
            self.live_state.finish(status="failed", summary=f"Hierarchical evaluation aborted: {exc}")
            raise
