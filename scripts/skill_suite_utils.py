#!/usr/bin/env python3
"""Helpers shared by skill-suite audit and benchmark tooling."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "registry"
TESTS = ROOT / "tests"
MAKEFILE = ROOT / "Makefile"


def load_jsonl(path: Path) -> list[dict]:
    records: list[dict] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        records.append(json.loads(line))
    return records


def load_skills() -> list[dict]:
    return load_jsonl(REGISTRY / "skills.jsonl")


def load_resources() -> dict[str, dict]:
    return {
        record["resource_id"]: record
        for record in load_jsonl(REGISTRY / "resources_dedup.jsonl")
    }


def load_metadata(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def count_asset_files(skill_path: Path) -> int:
    assets_dir = skill_path / "assets"
    if not assets_dir.exists():
        return 0
    return sum(
        1
        for path in assets_dir.rglob("*")
        if path.is_file() and path.name != ".gitkeep"
    )


def local_test_files(skill_path: Path) -> list[Path]:
    tests_dir = skill_path / "tests"
    if not tests_dir.exists():
        return []
    return sorted(tests_dir.glob("test_*.py"))


def parse_make_targets(path: Path = MAKEFILE) -> dict[str, list[str]]:
    targets: dict[str, list[str]] = {}
    current_target: str | None = None
    target_pattern = re.compile(r"^([A-Za-z0-9_.-]+):")
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip("\n")
        match = target_pattern.match(line)
        if match and not line.startswith("\t"):
            target = match.group(1)
            current_target = target if target.startswith("smoke-") else None
            if current_target is not None:
                targets.setdefault(current_target, [])
            continue
        if current_target is not None and raw_line.startswith("\t"):
            targets[current_target].append(raw_line.strip())
    return targets


def map_skill_to_smoke_targets(skills: list[dict], smoke_targets: dict[str, list[str]]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = defaultdict(list)
    for target, recipe_lines in smoke_targets.items():
        recipe_text = "\n".join(recipe_lines)
        for skill in skills:
            if skill["path"] in recipe_text or skill["slug"] in recipe_text:
                mapping[skill["skill_id"]].append(target)
    if "smoke-frontier-generated" in smoke_targets:
        for skill in skills:
            if "frontier-closure" in skill.get("tags", []):
                mapping[skill["skill_id"]].append("smoke-frontier-generated")
    return dict(mapping)


def build_repo_test_reference_map(skills: list[dict]) -> dict[str, list[str]]:
    references: dict[str, list[str]] = {skill["skill_id"]: [] for skill in skills}
    test_files = sorted(TESTS.rglob("test_*.py"))
    texts = {path: path.read_text(encoding="utf-8") for path in test_files}
    for skill in skills:
        skill_id = skill["skill_id"]
        for path, text in texts.items():
            if (
                skill["slug"] in text
                or skill["path"] in text
                or skill_id in text
                or ("frontier-closure" in skill.get("tags", []) and "frontier-closure" in text)
            ):
                references[skill_id].append(str(path.relative_to(ROOT)))
    return references
