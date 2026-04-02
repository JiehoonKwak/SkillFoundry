"""Prompt builders for codex-exec stages."""

from __future__ import annotations

import json

from .models import FocusLeaf, RepositorySummary

DEFAULT_CYCLE_STAGES = (
    "tree_check",
    "resource_search",
    "skill_build",
    "skill_test",
    "refresh",
)

SUPPORTED_STAGES = set(DEFAULT_CYCLE_STAGES)
PARALLEL_LEAF_STAGES = {"resource_search", "skill_build", "skill_test"}

STAGE_OBJECTIVES = {
    "tree_check": (
        "Audit the taxonomy, site tree, and registry coverage. Refine the focus list, update planning files, "
        "and identify the highest-value leaves for the next pass."
    ),
    "resource_search": (
        "Search broadly for canonical resources for the selected leaves, update the resource registries, "
        "deduplicate new entries, and record provenance in experiments.md and progress files."
    ),
    "skill_build": (
        "Implement or improve concrete skills for the selected leaves. Reuse existing scripts, registries, "
        "and skill patterns instead of creating parallel abstractions."
    ),
    "skill_test": (
        "Design and run repository-level smoke, regression, integration, and Slurm-facing checks for the updated skills. "
        "Only promote verification statuses after the corresponding checks pass."
    ),
    "refresh": (
        "Refresh site data, registry reports, and per-run reports so the repository remains internally consistent "
        "after the changes from earlier stages."
    ),
}

EVALUATION_REUSE_SURFACE = [
    "`scripts/validate_repository.py`",
    "`scripts/build_site.py`",
    "`scripts/audit_skill_suite.py`",
    "`scripts/run_skill_smoke_matrix.py`",
    "`scripts/benchmark_skill_advantage.py`",
    "`registry/skills.jsonl`",
    "`registry/resources_dedup.jsonl`",
    "`experiments.md`",
    "`README.md`",
    "`task_plan.md`, `findings.md`, `progress.md`",
]


def focus_leaf_payload(focus_leaves: list[FocusLeaf]) -> list[dict]:
    return [leaf.to_dict() for leaf in focus_leaves]


def stage_output_schema(stage: str) -> dict:
    if stage not in SUPPORTED_STAGES and stage != "design_skill":
        raise ValueError(f"Unsupported stage: {stage}")
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "stage": {"type": "string", "enum": [stage]},
            "summary": {"type": "string"},
            "selected_topics": {"type": "array", "items": {"type": "string"}},
            "skills_touched": {"type": "array", "items": {"type": "string"}},
            "resources_touched": {"type": "array", "items": {"type": "string"}},
            "tests_run": {"type": "array", "items": {"type": "string"}},
            "slurm_jobs": {"type": "array", "items": {"type": "string"}},
            "repo_changes": {"type": "array", "items": {"type": "string"}},
            "blockers": {"type": "array", "items": {"type": "string"}},
            "next_steps": {"type": "array", "items": {"type": "string"}},
        },
        "required": [
            "stage",
            "summary",
            "selected_topics",
            "skills_touched",
            "resources_touched",
            "tests_run",
            "slurm_jobs",
            "repo_changes",
            "blockers",
            "next_steps",
        ],
    }


def build_stage_prompt(
    *,
    stage: str,
    repo_summary: RepositorySummary,
    focus_leaves: list[FocusLeaf],
    extra_context: str | None = None,
) -> str:
    if stage not in SUPPORTED_STAGES:
        raise ValueError(f"Unsupported stage: {stage}")
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "Follow the workflow contract in `experiments.md` and reuse the existing repository scripts instead of building parallel registries.",
        "",
        f"STAGE: {stage}",
        "MODE: cycle",
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "FOCUS_LEAVES_JSON:",
        json.dumps(focus_leaf_payload(focus_leaves), indent=2, sort_keys=True),
        "",
        "STAGE_OBJECTIVE:",
        STAGE_OBJECTIVES[stage],
        "",
        "REQUIRED_REUSE_SURFACE:",
        "- `scripts/build_site.py`",
        "- `scripts/validate_repository.py`",
        "- `scripts/audit_skill_suite.py`",
        "- `scripts/run_skill_smoke_matrix.py`",
        "- `scripts/benchmark_skill_advantage.py`",
        "- `scripts/close_frontier_leaves.py`",
        "- `registry/taxonomy.yaml`",
        "- `registry/resources_dedup.jsonl`",
        "- `registry/skills.jsonl`",
        "- `site/tree.json`",
        "- `task_plan.md`, `findings.md`, `progress.md`",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Do the repository work for this stage directly in the current checkout when appropriate.",
        "- Keep verification labels honest.",
        "- Return only a JSON object that matches the provided schema.",
    ]
    return "\n".join(lines) + "\n"


def build_parallel_leaf_stage_prompt(
    *,
    stage: str,
    repo_summary: RepositorySummary,
    focus_leaf: FocusLeaf,
    artifact_dir: str,
    extra_context: str | None = None,
) -> str:
    if stage not in PARALLEL_LEAF_STAGES:
        raise ValueError(f"Unsupported parallel leaf stage: {stage}")
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside an isolated workspace copy of the SciSkillUniverse repository.",
        "This is a parallel worker stage for a single taxonomy leaf.",
        "Follow `experiments.md`, but stay scoped to the assigned leaf and avoid shared global files unless absolutely necessary.",
        "",
        f"STAGE: {stage}",
        "MODE: parallel_leaf_stage",
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "TARGET_FOCUS_LEAF_JSON:",
        json.dumps(focus_leaf.to_dict(), indent=2, sort_keys=True),
        "",
        "STAGE_OBJECTIVE:",
        STAGE_OBJECTIVES[stage],
        "",
        "ARTIFACT_DIRECTORY:",
        artifact_dir,
        "",
        "WORKER_RULES:",
        "- Treat this workspace as isolated. You may edit skill-local files and leaf-specific tests here.",
        "- Prefer changes under `skills/`, leaf-specific files under `tests/`, and leaf-specific reports inside the provided artifact directory.",
        "- Do not edit shared global files such as `registry/*.jsonl`, `registry/aliases.json`, `site/*.json`, `README.md`, `experiments.md`, or the planning files during this worker stage.",
        "- If shared-file follow-up is needed, record it in `repo_changes`, `blockers`, or `next_steps` instead of editing those files directly here.",
        "- Keep changes focused on this target leaf; do not broaden scope to unrelated domains.",
        "- When creating new file names, sanitize leaf labels into safe ASCII slugs: use dash-separated lowercase for skill directories and snake_case for `tests/test_*.py` files.",
        "- Avoid raw punctuation from taxonomy labels in filenames; convert characters such as `&`, `/`, or spaces into `-` or `_` first.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object that matches the provided schema.",
        "- Use `repo_changes` to list the concrete leaf-owned files you created or modified in this worker workspace.",
    ]
    return "\n".join(lines) + "\n"


def build_design_skill_prompt(
    *,
    task_prompt: str,
    repo_summary: RepositorySummary,
    focus_leaves: list[FocusLeaf],
    extra_context: str | None = None,
) -> str:
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "Follow the workflow contract in `experiments.md`.",
        "This mode is for a user-specified task: search resources first when the request needs additional materials, then implement and test the resulting skill.",
        "",
        "STAGE: design_skill",
        "MODE: design_skill",
        "",
        "USER_TASK_PROMPT:",
        task_prompt,
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "FOCUS_LEAVES_JSON:",
        json.dumps(focus_leaf_payload(focus_leaves), indent=2, sort_keys=True),
        "",
        "REQUIRED_BEHAVIOR:",
        "- Search for canonical papers, repositories, docs, notebooks, workflows, or package references if the task needs them.",
        "- Build or refine a single concrete skill path for the user task.",
        "- Add or update repository-level tests and, if relevant, a Slurm exercise path.",
        "- Update registries, reports, and `experiments.md` if repository state changes.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object that matches the provided schema.",
    ]
    return "\n".join(lines) + "\n"


def layer1_fix_output_schema() -> dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "stage": {"type": "string", "enum": ["layer1_fix"]},
            "skill_slug": {"type": "string"},
            "summary": {"type": "string"},
            "debug_findings": {"type": "array", "items": {"type": "string"}},
            "fixes_applied": {"type": "array", "items": {"type": "string"}},
            "tests_run": {"type": "array", "items": {"type": "string"}},
            "repo_changes": {"type": "array", "items": {"type": "string"}},
            "blockers": {"type": "array", "items": {"type": "string"}},
            "next_steps": {"type": "array", "items": {"type": "string"}},
        },
        "required": [
            "stage",
            "skill_slug",
            "summary",
            "debug_findings",
            "fixes_applied",
            "tests_run",
            "repo_changes",
            "blockers",
            "next_steps",
        ],
    }


def layer2_benchmark_output_schema() -> dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "stage": {"type": "string", "enum": ["layer2_benchmark"]},
            "skill_slug": {"type": "string"},
            "summary": {"type": "string"},
            "benchmark_cases": {"type": "array", "items": {"type": "string"}},
            "with_skill_score": {"type": "number"},
            "without_skill_score": {"type": "number"},
            "advantage_score": {"type": "number"},
            "meets_expectation": {"type": "boolean"},
            "optimization_targets": {"type": "array", "items": {"type": "string"}},
            "tests_run": {"type": "array", "items": {"type": "string"}},
            "repo_changes": {"type": "array", "items": {"type": "string"}},
            "blockers": {"type": "array", "items": {"type": "string"}},
            "next_steps": {"type": "array", "items": {"type": "string"}},
        },
        "required": [
            "stage",
            "skill_slug",
            "summary",
            "benchmark_cases",
            "with_skill_score",
            "without_skill_score",
            "advantage_score",
            "meets_expectation",
            "optimization_targets",
            "tests_run",
            "repo_changes",
            "blockers",
            "next_steps",
        ],
    }


def layer2_optimize_output_schema() -> dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "stage": {"type": "string", "enum": ["layer2_optimize"]},
            "skill_slug": {"type": "string"},
            "summary": {"type": "string"},
            "optimization_findings": {"type": "array", "items": {"type": "string"}},
            "optimizations_applied": {"type": "array", "items": {"type": "string"}},
            "tests_run": {"type": "array", "items": {"type": "string"}},
            "repo_changes": {"type": "array", "items": {"type": "string"}},
            "blockers": {"type": "array", "items": {"type": "string"}},
            "next_steps": {"type": "array", "items": {"type": "string"}},
        },
        "required": [
            "stage",
            "skill_slug",
            "summary",
            "optimization_findings",
            "optimizations_applied",
            "tests_run",
            "repo_changes",
            "blockers",
            "next_steps",
        ],
    }


def novelty_check_output_schema() -> dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "stage": {"type": "string", "enum": ["novelty_check"]},
            "skill_slug": {"type": "string"},
            "summary": {"type": "string"},
            "novelty_score": {"type": "number"},
            "local_overlap_slugs": {"type": "array", "items": {"type": "string"}},
            "external_overlap_titles": {"type": "array", "items": {"type": "string"}},
            "web_sources": {"type": "array", "items": {"type": "string"}},
            "needs_merge_review": {"type": "boolean"},
            "blockers": {"type": "array", "items": {"type": "string"}},
            "next_steps": {"type": "array", "items": {"type": "string"}},
        },
        "required": [
            "stage",
            "skill_slug",
            "summary",
            "novelty_score",
            "local_overlap_slugs",
            "external_overlap_titles",
            "web_sources",
            "needs_merge_review",
            "blockers",
            "next_steps",
        ],
    }


def build_layer1_fix_prompt(
    *,
    skill: dict,
    failure_record: dict,
    repo_summary: RepositorySummary,
    artifact_dir: str,
    extra_context: str | None = None,
) -> str:
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "This is hierarchical evaluation layer 1: correctness repair after a sandbox/Slurm-facing check failed.",
        "Diagnose the failure, repair the skill or its tests if needed, keep statuses honest, and return structured debug information.",
        "",
        "STAGE: layer1_fix",
        "MODE: evaluate_skills",
        "",
        "SKILL_JSON:",
        json.dumps(skill, indent=2, sort_keys=True),
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "FAILED_CHECK_JSON:",
        json.dumps(failure_record, indent=2, sort_keys=True),
        "",
        "ARTIFACT_DIR:",
        artifact_dir,
        "",
        "REQUIRED_REUSE_SURFACE:",
        *[f"- {item}" for item in EVALUATION_REUSE_SURFACE],
        "",
        "REQUIRED_BEHAVIOR:",
        "- Read the failure logs and inspect the skill implementation before changing code.",
        "- Apply a targeted fix instead of broad unrelated edits.",
        "- If the failure comes from missing tests or smoke coverage, repair the verification path too.",
        "- Write any temporary debug artifacts under the provided artifact directory.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object matching the schema.",
    ]
    return "\n".join(lines) + "\n"


def build_layer2_benchmark_prompt(
    *,
    skill: dict,
    repo_summary: RepositorySummary,
    artifact_dir: str,
    extra_context: str | None = None,
) -> str:
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "This is hierarchical evaluation layer 2: design and execute a task-specific benchmark that compares Codex with the skill versus without the skill.",
        "The goal is to prove that the maintained skill wrapper improves execution quality, output quality, or task completeness compared with an ad hoc baseline.",
        "",
        "STAGE: layer2_benchmark",
        "MODE: evaluate_skills",
        "",
        "SKILL_JSON:",
        json.dumps(skill, indent=2, sort_keys=True),
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "ARTIFACT_DIR:",
        artifact_dir,
        "",
        "REQUIRED_REUSE_SURFACE:",
        *[f"- {item}" for item in EVALUATION_REUSE_SURFACE],
        "",
        "REQUIRED_BEHAVIOR:",
        "- Design 1 to 3 concrete benchmark cases for this skill.",
        "- Execute a with-skill path and a no-skill baseline path for the cases when feasible.",
        "- Save benchmark notes or outputs under the provided artifact directory.",
        "- Use conservative scoring and only claim an advantage if the evidence is clear.",
        "- If the skill is not yet better than the baseline, identify optimization targets.",
        "",
        "SCORING_GUIDANCE:",
        "- `with_skill_score` and `without_skill_score` should be normalized to the range [0, 1].",
        "- `advantage_score` should equal `with_skill_score - without_skill_score`.",
        "- Set `meets_expectation` to true only when the maintained skill path is meaningfully better or more complete.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object matching the schema.",
    ]
    return "\n".join(lines) + "\n"


def build_layer2_optimize_prompt(
    *,
    skill: dict,
    benchmark_result: dict,
    repo_summary: RepositorySummary,
    artifact_dir: str,
    extra_context: str | None = None,
) -> str:
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "This is hierarchical evaluation layer 2 feedback: the skill did not outperform the no-skill baseline strongly enough.",
        "Improve the skill, its tests, or its promptable contract so the next benchmark run has a better chance of succeeding.",
        "",
        "STAGE: layer2_optimize",
        "MODE: evaluate_skills",
        "",
        "SKILL_JSON:",
        json.dumps(skill, indent=2, sort_keys=True),
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "BENCHMARK_RESULT_JSON:",
        json.dumps(benchmark_result, indent=2, sort_keys=True),
        "",
        "ARTIFACT_DIR:",
        artifact_dir,
        "",
        "REQUIRED_REUSE_SURFACE:",
        *[f"- {item}" for item in EVALUATION_REUSE_SURFACE],
        "",
        "REQUIRED_BEHAVIOR:",
        "- Focus on the specific deficits reported by the benchmark.",
        "- Preserve or improve correctness while optimizing for benchmark advantage.",
        "- Update tests or examples if they are the limiting factor.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object matching the schema.",
    ]
    return "\n".join(lines) + "\n"


def build_novelty_check_prompt(
    *,
    skill: dict,
    local_candidates: list[dict],
    repo_summary: RepositorySummary,
    artifact_dir: str,
    extra_context: str | None = None,
) -> str:
    context_block = extra_context.strip() if extra_context else "None."
    lines = [
        "You are operating inside the SciSkillUniverse repository.",
        "Run a novelty and overlap check for the selected skill.",
        "Use both the local registry and web search to check overlap with existing resources and prior skill libraries.",
        "",
        "STAGE: novelty_check",
        "MODE: evaluate_skills",
        "",
        "SKILL_JSON:",
        json.dumps(skill, indent=2, sort_keys=True),
        "",
        "CURRENT_REPOSITORY_SUMMARY_JSON:",
        json.dumps(repo_summary.to_dict(), indent=2, sort_keys=True),
        "",
        "LOCAL_SIMILARITY_CANDIDATES_JSON:",
        json.dumps(local_candidates, indent=2, sort_keys=True),
        "",
        "ARTIFACT_DIR:",
        artifact_dir,
        "",
        "REQUIRED_BEHAVIOR:",
        "- Check overlap with similar local skills first.",
        "- Use web search to inspect external overlapping skill or tool ecosystems when relevant.",
        "- Explicitly check overlap against ToolUniverse: https://github.com/mims-harvard/ToolUniverse",
        "- Save any compact notes under the provided artifact directory if useful.",
        "- Recommend merge review when the skill appears duplicative or insufficiently differentiated.",
        "",
        "SCORING_GUIDANCE:",
        "- `novelty_score` should be in the range [0, 5], where 0 means nearly duplicate and 5 means clearly differentiated.",
        "",
        "ADDITIONAL_CONTEXT:",
        context_block,
        "",
        "RESPONSE_CONTRACT:",
        "- Return only a JSON object matching the schema.",
    ]
    return "\n".join(lines) + "\n"
