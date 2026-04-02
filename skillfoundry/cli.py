"""CLI for the codex-exec automation framework."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .campaign import FrameworkCampaignRunner
from .orchestrator import SkillAutomationFramework
from .prompts import DEFAULT_CYCLE_STAGES


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--codex-bin", default="codex")
    parser.add_argument("--state-dir", type=Path, default=None)
    parser.add_argument("--reports-dir", type=Path, default=None)
    parser.add_argument("--workspace-root", type=Path, default=None)
    parser.add_argument("--model", default=None)
    parser.add_argument("--reasoning-effort", default=None)
    parser.add_argument("--profile", default=None)
    parser.add_argument("--full-auto", action="store_true")
    parser.add_argument(
        "--codex-config-override",
        action="append",
        default=[],
        help="Repeatable raw `codex exec -c key=value` override, e.g. --codex-config-override 'features.apps=false'",
    )
    parser.add_argument("--stage-timeout", type=int, default=3600)
    parser.add_argument("--verification-timeout", type=int, default=3600)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    status = subparsers.add_parser("status", help="Show the current repository summary and focus leaves.")
    status.add_argument("--focus-limit", type=int, default=8)
    status.add_argument("--focus-term", action="append", default=[], help="Limit focus leaves to matching domains / topics.")

    cycle = subparsers.add_parser("cycle", help="Run the staged codex-exec automation loop.")
    cycle.add_argument("--loops", type=int, default=1)
    cycle.add_argument("--focus-limit", type=int, default=8)
    cycle.add_argument("--focus-term", action="append", default=[], help="Limit focus leaves to matching domains / topics.")
    cycle.add_argument("--stage-workers", type=int, default=1, help="Concurrent codex-exec workers for leaf-parallel stages.")
    cycle.add_argument("--background-validation-limit", type=int, default=0, help="Run background skill validation for this many skills while worker stages execute.")
    cycle.add_argument("--background-validation-workers", type=int, default=1, help="Concurrent workers for the background validation lane.")
    cycle.add_argument("--verification-mode", default="standard", choices=["none", "validate", "standard", "full", "audit"])
    cycle.add_argument("--label", default="cycle")
    cycle.add_argument("--extra-context", default=None)
    cycle.add_argument(
        "--stages",
        default=",".join(DEFAULT_CYCLE_STAGES),
        help="Comma-separated stage list.",
    )

    design = subparsers.add_parser("design-skill", help="Design a skill for a specific task prompt.")
    design.add_argument("--prompt", required=True)
    design.add_argument("--focus-limit", type=int, default=8)
    design.add_argument("--focus-term", action="append", default=[], help="Limit focus leaves to matching domains / topics.")
    design.add_argument("--verification-mode", default="standard", choices=["none", "validate", "standard", "full", "audit"])
    design.add_argument("--label", default="design-skill")
    design.add_argument("--extra-context", default=None)

    evaluate = subparsers.add_parser(
        "evaluate-skills",
        help="Run hierarchical correctness, benchmark, and novelty evaluation for existing skills.",
    )
    evaluate.add_argument("--skill-slug", action="append", default=[])
    evaluate.add_argument("--all", action="store_true", help="Evaluate the full registered skill library.")
    evaluate.add_argument("--limit", type=int, default=8, help="Number of skills to evaluate when --all is not set.")
    evaluate.add_argument("--verification-mode", default="standard", choices=["none", "validate", "standard", "full", "audit"])
    evaluate.add_argument("--label", default="evaluate-skills")
    evaluate.add_argument("--workers", type=int, default=1, help="Concurrent workers per evaluation layer.")
    evaluate.add_argument("--layer1-fix-attempts", type=int, default=1)
    evaluate.add_argument("--layer2-optimize-attempts", type=int, default=1)
    evaluate.add_argument("--skip-novelty-check", action="store_true")
    evaluate.add_argument("--extra-context", default=None)

    campaign = subparsers.add_parser("campaign", help="Run a checkpointable long campaign across targeted domains.")
    campaign.add_argument("--label", default="domain-campaign")
    campaign.add_argument("--focus-term", action="append", default=[], required=True, help="Target domain/topic terms for discovery focus.")
    campaign.add_argument("--focus-limit", type=int, default=12)
    campaign.add_argument("--stage-workers", type=int, default=6)
    campaign.add_argument("--background-validation-limit", type=int, default=24)
    campaign.add_argument("--background-validation-workers", type=int, default=6)
    campaign.add_argument("--evaluation-batch-size", type=int, default=24)
    campaign.add_argument("--evaluation-workers", type=int, default=6)
    campaign.add_argument("--layer1-fix-attempts", type=int, default=1)
    campaign.add_argument("--layer2-optimize-attempts", type=int, default=1)
    campaign.add_argument("--verification-mode", default="standard", choices=["none", "validate", "standard", "full", "audit"])
    campaign.add_argument("--full-verification-every", type=int, default=0)
    campaign.add_argument("--max-iterations", type=int, default=1000)
    campaign.add_argument("--max-runtime-minutes", type=int, default=450)
    campaign.add_argument("--stop-buffer-minutes", type=int, default=20)
    campaign.add_argument("--extra-context", default=None)

    campaign_status = subparsers.add_parser("campaign-status", help="Print the latest checkpoint summary for a campaign.")
    campaign_status.add_argument("--label", required=True)

    return parser


def framework_from_args(args: argparse.Namespace) -> SkillAutomationFramework:
    return SkillAutomationFramework(
        args.repo_root,
        codex_bin=args.codex_bin,
        state_dir=args.state_dir,
        reports_dir=args.reports_dir,
        workspace_root=args.workspace_root,
        model=args.model,
        reasoning_effort=args.reasoning_effort,
        profile=args.profile,
        full_auto=args.full_auto,
        codex_config_overrides=args.codex_config_override,
        stage_timeout=args.stage_timeout,
        verification_timeout=args.verification_timeout,
    )


def summary_from_payload(payload: dict) -> dict | None:
    if "evaluation_summary" in payload and isinstance(payload["evaluation_summary"], dict):
        return payload["evaluation_summary"]
    if "summary" in payload and isinstance(payload["summary"], dict):
        return payload["summary"]
    if "campaign_summary" in payload and isinstance(payload["campaign_summary"], dict):
        return payload["campaign_summary"]
    if "final_status" in payload and isinstance(payload["final_status"], dict):
        final_status = payload["final_status"]
        if "summary" in final_status:
            return final_status["summary"]
        return final_status
    if "initial_status" in payload and isinstance(payload["initial_status"], dict):
        initial_status = payload["initial_status"]
        if "summary" in initial_status:
            return initial_status["summary"]
        return initial_status
    return None


def print_payload(payload: dict, emit_json: bool) -> None:
    if emit_json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    summary = summary_from_payload(payload)
    if summary:
        if {"resource_count", "skill_count", "covered_leaf_count", "frontier_leaf_count", "todo_leaf_count"} <= set(summary):
            print(
                "Resources={resource_count} Skills={skill_count} CoveredLeaves={covered_leaf_count} "
                "FrontierLeaves={frontier_leaf_count} TodoLeaves={todo_leaf_count}".format(**summary)
            )
        elif {"selected_skill_count", "layer1_passed", "layer2_advantage_passed"} <= set(summary):
            print(
                "SelectedSkills={selected_skill_count} "
                "Layer1Passed={layer1_passed} Layer1Failed={layer1_failed} "
                "Layer2Passed={layer2_advantage_passed} Layer2Failed={layer2_advantage_failed} "
                "NoveltyMergeReview={novelty_merge_review_count}".format(**summary)
            )
        elif {"designed_skill_count", "new_skill_count", "layer1_pass_count", "layer2_pass_count"} <= set(summary):
            print(
                "CampaignStatus={status} Iterations={iterations_completed} "
                "Designed={designed_skill_count} New={new_skill_count} "
                "Layer1Passed={layer1_pass_count} Layer2Passed={layer2_pass_count} "
                "PendingEval={pending_evaluation_count}".format(**summary)
            )
    focus_leaves = payload.get("focus_leaves")
    if focus_leaves:
        print("Focus leaves:")
        for leaf in focus_leaves:
            print(
                f"- {leaf['domain_slug']} / {leaf['topic_slug']}: "
                f"{leaf['coverage_status']} (resources={leaf['resource_count']}, skills={leaf['skill_count']}, verified={leaf['verified_skill_count']})"
            )
    elif "run_dir" in payload:
        print(f"Run artifacts: {payload['run_dir']}")
        print(f"Verification records: {len(payload.get('verification', []))}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    framework = framework_from_args(args)
    campaign_runner = FrameworkCampaignRunner(framework)

    if args.command == "status":
        payload = framework.status(focus_limit=args.focus_limit, focus_terms=args.focus_term or None)
    elif args.command == "cycle":
        stages = [stage.strip() for stage in args.stages.split(",") if stage.strip()]
        payload = framework.run_cycle(
            loops=args.loops,
            stages=stages,
            focus_limit=args.focus_limit,
            focus_terms=args.focus_term or None,
            stage_workers=args.stage_workers,
            background_validation_limit=args.background_validation_limit,
            background_validation_workers=args.background_validation_workers,
            verification_mode=args.verification_mode,
            run_label=args.label,
            extra_context=args.extra_context,
        )
    elif args.command == "evaluate-skills":
        payload = framework.evaluate_skills(
            skill_slugs=args.skill_slug or None,
            include_all=args.all,
            limit=args.limit,
            workers=args.workers,
            verification_mode=args.verification_mode,
            run_label=args.label,
            layer1_fix_attempts=args.layer1_fix_attempts,
            layer2_optimize_attempts=args.layer2_optimize_attempts,
            novelty_check=not args.skip_novelty_check,
            extra_context=args.extra_context,
        )
    elif args.command == "campaign":
        payload = campaign_runner.run_campaign(
            label=args.label,
            focus_terms=args.focus_term,
            cycle_focus_limit=args.focus_limit,
            cycle_stage_workers=args.stage_workers,
            background_validation_limit=args.background_validation_limit,
            background_validation_workers=args.background_validation_workers,
            evaluation_batch_size=args.evaluation_batch_size,
            evaluation_workers=args.evaluation_workers,
            layer1_fix_attempts=args.layer1_fix_attempts,
            layer2_optimize_attempts=args.layer2_optimize_attempts,
            verification_mode=args.verification_mode,
            full_verification_every=args.full_verification_every,
            max_iterations=args.max_iterations,
            max_runtime_minutes=args.max_runtime_minutes,
            stop_buffer_minutes=args.stop_buffer_minutes,
            extra_context=args.extra_context,
        )
    elif args.command == "campaign-status":
        payload = campaign_runner.campaign_status(label=args.label)
    else:
        payload = framework.design_skill(
            task_prompt=args.prompt,
            focus_limit=args.focus_limit,
            focus_terms=args.focus_term or None,
            verification_mode=args.verification_mode,
            run_label=args.label,
            extra_context=args.extra_context,
        )

    print_payload(payload, args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
