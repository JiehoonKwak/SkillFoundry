#!/usr/bin/env python3
"""Audit structure, provenance, and coverage across the registered skill suite."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from skill_suite_utils import (
    ROOT,
    build_repo_test_reference_map,
    count_asset_files,
    load_metadata,
    load_resources,
    load_skills,
    local_test_files,
    map_skill_to_smoke_targets,
    parse_make_targets,
)


REQUIRED_NAMES = ["SKILL.md", "metadata.yaml", "refs.md", "scripts", "examples", "tests", "assets"]


def audit_skill(skill: dict, resources: dict[str, dict], smoke_map: dict[str, list[str]], repo_test_refs: dict[str, list[str]]) -> dict:
    skill_path = ROOT / skill["path"]
    local_tests = local_test_files(skill_path)
    repo_tests = repo_test_refs.get(skill["skill_id"], [])
    smoke_targets = smoke_map.get(skill["skill_id"], [])
    required_paths = {name: (skill_path / name).exists() for name in REQUIRED_NAMES}
    hard_failures: list[str] = []
    warnings: list[str] = []

    if not skill_path.exists():
        hard_failures.append(f"missing skill path: {skill['path']}")
        metadata = {}
    else:
        metadata = load_metadata(skill_path / "metadata.yaml")
        if metadata.get("slug") != skill["slug"]:
            hard_failures.append("metadata slug mismatch")
        if metadata.get("name") != skill["name"]:
            hard_failures.append("metadata name mismatch")
        if metadata.get("status") != skill["status"]:
            hard_failures.append("metadata status mismatch")

    for name, exists in required_paths.items():
        if not exists:
            hard_failures.append(f"missing required path: {name}")

    unknown_resources = [resource_id for resource_id in skill["source_resource_ids"] if resource_id not in resources]
    if unknown_resources:
        hard_failures.append(f"unknown source resources: {', '.join(unknown_resources)}")
    if not skill.get("source_resource_ids"):
        hard_failures.append("no source resources declared")
    if not skill.get("test_commands"):
        hard_failures.append("no test_commands declared")
    if not smoke_targets:
        hard_failures.append("no Makefile smoke target")
    if not local_tests:
        warnings.append("no skill-local tests")
    if not repo_tests:
        warnings.append("no repository-level test references")
    if count_asset_files(skill_path) == 0:
        warnings.append("no asset files")

    return {
        "skill_id": skill["skill_id"],
        "slug": skill["slug"],
        "path": skill["path"],
        "status": skill["status"],
        "skill_type": skill["skill_type"],
        "domain": skill["topic_path"][0],
        "source_resource_count": len(skill["source_resource_ids"]),
        "asset_file_count": count_asset_files(skill_path),
        "local_test_files": [str(path.relative_to(ROOT)) for path in local_tests],
        "local_test_count": len(local_tests),
        "repo_test_refs": repo_tests,
        "repo_test_ref_count": len(repo_tests),
        "smoke_targets": smoke_targets,
        "required_paths": required_paths,
        "hard_failures": hard_failures,
        "warnings": warnings,
    }


def build_summary(records: list[dict], skills: list[dict]) -> dict:
    no_local_tests = [record["skill_id"] for record in records if "no skill-local tests" in record["warnings"]]
    summary = {
        "skill_count": len(records),
        "status_counts": dict(Counter(skill["status"] for skill in skills)),
        "domain_counts": dict(Counter(skill["topic_path"][0] for skill in skills)),
        "skill_type_counts": dict(Counter(skill["skill_type"] for skill in skills)),
        "skills_with_smoke_target": sum(1 for record in records if record["smoke_targets"]),
        "skills_with_local_tests": sum(1 for record in records if record["local_test_count"] > 0),
        "skills_with_repo_test_refs": sum(1 for record in records if record["repo_test_ref_count"] > 0),
        "skills_with_any_test_refs": sum(
            1
            for record in records
            if record["local_test_count"] > 0 or record["repo_test_ref_count"] > 0
        ),
        "skills_with_assets": sum(1 for record in records if record["asset_file_count"] > 0),
        "hard_failure_count": sum(len(record["hard_failures"]) for record in records),
        "warning_count": sum(len(record["warnings"]) for record in records),
        "skills_without_local_tests": no_local_tests,
        "skills_with_hard_failures": [
            record["skill_id"] for record in records if record["hard_failures"]
        ],
    }
    return summary


def write_json(data: dict, path: Path | None) -> None:
    text = json.dumps(data, indent=2, sort_keys=True)
    if path is None:
        print(text)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n", encoding="utf-8")


def write_markdown(summary: dict, records: list[dict], path: Path | None) -> None:
    if path is None:
        return
    lines = [
        "# Skill Suite Audit",
        "",
        f"- Skills audited: `{summary['skill_count']}`",
        f"- Skills with Makefile smoke coverage: `{summary['skills_with_smoke_target']}`",
        f"- Skills with local tests: `{summary['skills_with_local_tests']}`",
        f"- Skills with repository-level test references: `{summary['skills_with_repo_test_refs']}`",
        f"- Skills with assets: `{summary['skills_with_assets']}`",
        f"- Hard failures: `{summary['hard_failure_count']}`",
        f"- Warnings: `{summary['warning_count']}`",
        "",
        "## Findings",
        "",
    ]
    if summary["hard_failure_count"] == 0:
        lines.append("- No hard structural failures were detected.")
    else:
        for record in records:
            for failure in record["hard_failures"]:
                lines.append(f"- `{record['skill_id']}`: {failure}")
    if summary["skills_without_local_tests"]:
        lines.append(
            f"- `{len(summary['skills_without_local_tests'])}` legacy skills still rely on repository-level coverage instead of skill-local tests."
        )
    lines.extend(
        [
            "",
            "## Skills Without Local Tests",
            "",
        ]
    )
    if summary["skills_without_local_tests"]:
        for skill_id in summary["skills_without_local_tests"]:
            lines.append(f"- `{skill_id}`")
    else:
        lines.append("- None")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output path.")
    parser.add_argument("--markdown-out", type=Path, default=None, help="Optional Markdown summary path.")
    args = parser.parse_args()

    skills = load_skills()
    resources = load_resources()
    smoke_targets = parse_make_targets()
    smoke_map = map_skill_to_smoke_targets(skills, smoke_targets)
    repo_test_refs = build_repo_test_reference_map(skills)

    records = [
        audit_skill(skill=skill, resources=resources, smoke_map=smoke_map, repo_test_refs=repo_test_refs)
        for skill in skills
    ]
    summary = build_summary(records, skills)
    payload = {"summary": summary, "skills": records}

    write_json(payload, args.json_out)
    write_markdown(summary, records, args.markdown_out)

    if summary["hard_failure_count"] > 0:
        print("Skill suite audit found hard failures.", file=sys.stderr)
        return 1
    print(
        "Skill suite audit passed: "
        f"{summary['skill_count']} skills, "
        f"{summary['skills_with_smoke_target']} smoke-covered, "
        f"{summary['skills_with_local_tests']} with local tests."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
