"""Repository access helpers for the automation framework."""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import time
from difflib import SequenceMatcher
from pathlib import Path
from types import ModuleType

from .models import FocusLeaf, RepositorySummary

VERIFIED_STATUSES = {"sandbox_verified", "slurm_verified"}


def load_script_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class SciSkillRepository:
    """Thin wrapper around the existing registry and helper scripts."""

    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.registry = self.root / "registry"
        self.site = self.root / "site"
        self.scratch = self.root / "scratch"
        self.reports = self.root / "reports"
        self._build_site = load_script_module(
            "sciskill_framework_build_site",
            self.root / "scripts" / "build_site.py",
        )
        self._skill_suite_utils = load_script_module(
            "sciskill_framework_skill_suite_utils",
            self.root / "scripts" / "skill_suite_utils.py",
        )

    def load_json(self, path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))

    def load_jsonl(self, path: Path) -> list[dict]:
        return self._build_site.load_jsonl(path)

    def load_skills(self) -> list[dict]:
        return self.load_jsonl(self.registry / "skills.jsonl")

    def load_resources(self) -> list[dict]:
        return self.load_jsonl(self.registry / "resources_dedup.jsonl")

    def load_taxonomy(self) -> dict:
        return self._build_site.load_json(self.registry / "taxonomy.yaml")

    def _normalize_focus_terms(self, focus_terms: list[str] | None) -> list[str]:
        return [term.strip().lower() for term in (focus_terms or []) if term and term.strip()]

    def _matches_focus_terms(self, values: list[object], focus_terms: list[str] | None) -> bool:
        normalized_terms = self._normalize_focus_terms(focus_terms)
        if not normalized_terms:
            return True
        raw_text = " ".join(str(value or "") for value in values).lower()
        slug_text = re.sub(r"[^a-z0-9]+", "-", raw_text).strip("-")
        plain_text = re.sub(r"[^a-z0-9]+", " ", raw_text).strip()
        for term in normalized_terms:
            term_slug = re.sub(r"[^a-z0-9]+", "-", term).strip("-")
            term_plain = re.sub(r"[^a-z0-9]+", " ", term).strip()
            if (term and term in raw_text) or (term_slug and term_slug in slug_text) or (term_plain and term_plain in plain_text):
                return True
        return False

    def skill_index(self) -> dict[str, dict]:
        return {skill["skill_id"]: skill for skill in self.load_skills()}

    def smoke_targets(self) -> dict[str, list[str]]:
        return self._skill_suite_utils.parse_make_targets(self.root / "Makefile")

    def smoke_map(self, skills: list[dict] | None = None) -> dict[str, list[str]]:
        records = skills if skills is not None else self.load_skills()
        return self._skill_suite_utils.map_skill_to_smoke_targets(records, self.smoke_targets())

    def choose_skill_target(self, skill: dict, smoke_map: dict[str, list[str]] | None = None) -> dict | None:
        resolved_map = smoke_map if smoke_map is not None else self.smoke_map([skill])
        smoke_targets = sorted(resolved_map.get(skill["skill_id"], []))
        if smoke_targets:
            target = smoke_targets[0]
            return {
                "kind": "make-target",
                "label": target,
                "command": ["make", target],
            }
        test_commands = skill.get("test_commands", [])
        if test_commands:
            command_text = test_commands[0]
            return {
                "kind": "registry-test-command",
                "label": command_text,
                "command": ["bash", "-lc", command_text],
            }
        return None

    def execute_skill_check(
        self,
        skill: dict,
        *,
        artifact_dir: Path | None = None,
        timeout: int = 1800,
        smoke_map: dict[str, list[str]] | None = None,
        progress_callback=None,
    ) -> dict:
        target = self.choose_skill_target(skill, smoke_map=smoke_map)
        if target is None:
            record = {
                "skill_id": skill["skill_id"],
                "slug": skill["slug"],
                "status": skill["status"],
                "target_kind": None,
                "target_label": None,
                "command": None,
                "returncode": 1,
                "duration_seconds": 0.0,
                "stdout_path": None,
                "stderr_path": None,
                "stdout_tail": [],
                "stderr_tail": ["No smoke target or registry test command could be resolved for this skill."],
            }
            if artifact_dir is not None:
                write_json(artifact_dir / "result.json", record)
            return record

        if artifact_dir is not None:
            artifact_dir.mkdir(parents=True, exist_ok=True)
        if progress_callback is not None:
            progress_callback(
                "start",
                {
                    "skill_id": skill["skill_id"],
                    "slug": skill["slug"],
                    "target_kind": target["kind"],
                    "target_label": target["label"],
                    "command": target["command"],
                },
            )

        started = time.monotonic()
        try:
            completed = subprocess.run(
                target["command"],
                cwd=self.root,
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
            duration = round(time.monotonic() - started, 3)
            stdout_text = completed.stdout
            stderr_text = completed.stderr
            returncode = completed.returncode
        except subprocess.TimeoutExpired as exc:
            duration = round(time.monotonic() - started, 3)
            stdout_text = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
            stderr_text = (exc.stderr or "") if isinstance(exc.stderr, str) else ""
            returncode = 124
        stdout_path = None
        stderr_path = None
        if artifact_dir is not None:
            stdout_path = artifact_dir / "stdout.txt"
            stderr_path = artifact_dir / "stderr.txt"
            stdout_path.write_text(stdout_text, encoding="utf-8")
            stderr_path.write_text(stderr_text, encoding="utf-8")

        record = {
            "skill_id": skill["skill_id"],
            "slug": skill["slug"],
            "status": skill["status"],
            "target_kind": target["kind"],
            "target_label": target["label"],
            "command": target["command"],
            "returncode": returncode,
            "duration_seconds": duration,
            "stdout_path": str(stdout_path) if stdout_path is not None else None,
            "stderr_path": str(stderr_path) if stderr_path is not None else None,
            "stdout_tail": stdout_text.strip().splitlines()[-20:],
            "stderr_tail": stderr_text.strip().splitlines()[-20:],
        }
        if artifact_dir is not None:
            write_json(artifact_dir / "result.json", record)
        if progress_callback is not None:
            progress_callback("finish", record)
        return record

    def select_skills(
        self,
        *,
        skill_slugs: list[str] | None = None,
        include_all: bool = False,
        limit: int = 8,
        focus_terms: list[str] | None = None,
    ) -> list[dict]:
        skills = self.load_skills()
        index = {
            value: skill
            for skill in skills
            for value in {skill["slug"], skill["skill_id"]}
        }
        if skill_slugs:
            missing = [slug for slug in skill_slugs if slug not in index]
            if missing:
                raise ValueError(f"Unknown skills requested: {', '.join(sorted(missing))}")
            return [index[slug] for slug in skill_slugs]
        if focus_terms:
            skills = [
                skill
                for skill in skills
                if self._matches_focus_terms(
                    [
                        skill.get("domain"),
                        skill.get("slug"),
                        skill.get("name"),
                        skill.get("description"),
                        *skill.get("tags", []),
                        *skill.get("topic_path", []),
                    ],
                    focus_terms,
                )
            ]
        if include_all:
            return sorted(skills, key=lambda skill: (skill["domain"], skill["slug"]))

        status_rank = {
            "partially_broken": 0,
            "stale": 1,
            "implemented": 2,
            "draft": 3,
            "idea": 4,
            "sandbox_verified": 5,
            "slurm_verified": 6,
            "deprecated": 7,
            "merged": 8,
        }
        ranked = sorted(
            skills,
            key=lambda skill: (
                status_rank.get(skill["status"], 99),
                skill["domain"],
                skill["slug"],
            ),
        )
        return ranked[:limit]

    def local_similarity_candidates(self, skill: dict, *, limit: int = 5) -> list[dict]:
        def slug_similarity(left: str, right: str) -> float:
            return SequenceMatcher(a=left, b=right).ratio()

        def tokenize(*values: object) -> set[str]:
            text = " ".join(str(value or "") for value in values).lower()
            return {token for token in re.findall(r"[a-z0-9]+", text) if len(token) > 2}

        topic_path = tuple(skill.get("topic_path", []))
        skill_tokens = tokenize(skill.get("slug"), skill.get("name"), skill.get("description"), *skill.get("tags", []))
        records: list[dict] = []
        for candidate in self.load_skills():
            if candidate["skill_id"] == skill["skill_id"]:
                continue
            candidate_tokens = tokenize(
                candidate.get("slug"),
                candidate.get("name"),
                candidate.get("description"),
                *candidate.get("tags", []),
            )
            shared_tokens = sorted(skill_tokens & candidate_tokens)
            same_domain = candidate.get("domain") == skill.get("domain")
            shared_topic = sorted(set(candidate.get("topic_path", [])) & set(topic_path))
            shared_resources = sorted(set(candidate.get("source_resource_ids", [])) & set(skill.get("source_resource_ids", [])))
            score = 0.0
            score += slug_similarity(skill["slug"], candidate["slug"]) * 2.5
            score += min(len(shared_tokens), 6) * 0.4
            score += len(shared_resources) * 1.25
            if same_domain:
                score += 1.0
            if shared_topic:
                score += 0.75 * len(shared_topic)
            if skill["slug"] in candidate.get("related_skills", []) or candidate["slug"] in skill.get("related_skills", []):
                score += 1.0
            if score <= 0:
                continue
            reasons = []
            if same_domain:
                reasons.append("same_domain")
            if shared_topic:
                reasons.append(f"shared_topic:{','.join(shared_topic)}")
            if shared_resources:
                reasons.append(f"shared_resources:{','.join(shared_resources)}")
            if shared_tokens:
                reasons.append(f"shared_terms:{','.join(shared_tokens[:6])}")
            records.append(
                {
                    "skill_id": candidate["skill_id"],
                    "slug": candidate["slug"],
                    "domain": candidate["domain"],
                    "status": candidate["status"],
                    "path": candidate["path"],
                    "score": round(score, 3),
                    "reasons": reasons,
                }
            )
        records.sort(key=lambda record: (-record["score"], record["domain"], record["slug"]))
        return records[:limit]

    def build_tree(self) -> dict:
        return self._build_site.build_tree(
            self.load_skills(),
            self.load_resources(),
            self.load_taxonomy(),
        )

    def ranked_focus_leaves(self, limit: int = 8, focus_terms: list[str] | None = None) -> list[FocusLeaf]:
        tree = self.build_tree()
        leaves: list[FocusLeaf] = []
        for domain in tree["children"]:
            domain_slug = domain["registry_domains"][0]
            for leaf in domain["children"]:
                skill_statuses = [skill["status"] for skill in leaf.get("skills", [])]
                verified_skill_count = sum(1 for status in skill_statuses if status in VERIFIED_STATUSES)
                leaves.append(
                    FocusLeaf(
                        taxonomy_key=domain["taxonomy_key"],
                        domain_slug=domain_slug,
                        domain_name=domain["name"],
                        leaf_name=leaf["name"],
                        topic_slug=leaf["topic_slug"],
                        coverage_status=leaf["coverage_status"],
                        skill_count=leaf["skill_count"],
                        resource_count=leaf["resource_count"],
                        verified_skill_count=verified_skill_count,
                        skill_slugs=[skill["slug"] for skill in leaf.get("skills", [])],
                        resource_ids=[resource["resource_id"] for resource in leaf.get("resources", [])],
                    )
                )
        if focus_terms:
            leaves = [
                leaf
                for leaf in leaves
                if self._matches_focus_terms(
                    [
                        leaf.taxonomy_key,
                        leaf.domain_slug,
                        leaf.domain_name,
                        leaf.leaf_name,
                        leaf.topic_slug,
                        *leaf.skill_slugs,
                        *leaf.resource_ids,
                    ],
                    focus_terms,
                )
            ]
        coverage_rank = {"todo": 0, "frontier": 1, "covered": 2}
        leaves.sort(
            key=lambda leaf: (
                coverage_rank.get(leaf.coverage_status, 3),
                leaf.verified_skill_count,
                -leaf.resource_count,
                leaf.skill_count,
                leaf.domain_slug,
                leaf.topic_slug,
            )
        )
        return leaves[:limit]

    def summary(self, focus_limit: int = 8, focus_terms: list[str] | None = None) -> RepositorySummary:
        tree = self.build_tree()
        covered_domain_count = sum(1 for node in tree["children"] if node["coverage_status"] == "covered")
        empty_domain_count = sum(1 for node in tree["children"] if node["coverage_status"] == "empty")
        return RepositorySummary(
            resource_count=tree["resource_count"],
            skill_count=tree["skill_count"],
            taxonomy_domain_count=tree["taxonomy_domain_count"],
            covered_leaf_count=tree["covered_leaf_count"],
            frontier_leaf_count=tree["frontier_leaf_count"],
            todo_leaf_count=tree["todo_leaf_count"],
            covered_domain_count=covered_domain_count,
            empty_domain_count=empty_domain_count,
            focus_leaves=self.ranked_focus_leaves(limit=focus_limit, focus_terms=focus_terms),
        )

    def verification_commands(self, mode: str) -> list[tuple[str, list[str]]]:
        commands = {
            "none": [],
            "validate": [("validate", ["python3", "scripts/validate_repository.py"])],
            "standard": [
                ("validate", ["python3", "scripts/validate_repository.py"]),
                ("build-site", ["python3", "scripts/build_site.py"]),
            ],
            "full": [
                ("validate", ["python3", "scripts/validate_repository.py"]),
                ("build-site", ["python3", "scripts/build_site.py"]),
                ("test", ["python3", "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]),
            ],
            "audit": [
                ("validate", ["python3", "scripts/validate_repository.py"]),
                ("build-site", ["python3", "scripts/build_site.py"]),
                ("test", ["python3", "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]),
                ("audit-skills", ["python3", "scripts/audit_skill_suite.py"]),
                ("smoke-matrix-dry-run", ["python3", "scripts/run_skill_smoke_matrix.py", "--dry-run"]),
            ],
        }
        if mode not in commands:
            raise ValueError(f"Unsupported verification mode: {mode}")
        return commands[mode]

    def execute_verification(
        self,
        mode: str,
        run_dir: Path | None = None,
        timeout: int = 3600,
        progress_callback=None,
    ) -> list[dict]:
        records: list[dict] = []
        log_dir = None if run_dir is None else run_dir / "verification"
        if log_dir is not None:
            log_dir.mkdir(parents=True, exist_ok=True)
        for label, command in self.verification_commands(mode):
            if progress_callback is not None:
                progress_callback(
                    "start",
                    {
                        "label": label,
                        "command": command,
                    },
                )
            started = time.monotonic()
            try:
                completed = subprocess.run(
                    command,
                    cwd=self.root,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )
                duration = round(time.monotonic() - started, 3)
                stdout_text = completed.stdout
                stderr_text = completed.stderr
                returncode = completed.returncode
            except subprocess.TimeoutExpired as exc:
                duration = round(time.monotonic() - started, 3)
                stdout_text = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
                stderr_text = (exc.stderr or "") if isinstance(exc.stderr, str) else ""
                returncode = 124
            stdout_path = None
            stderr_path = None
            if log_dir is not None:
                stdout_path = log_dir / f"{label}.stdout.txt"
                stderr_path = log_dir / f"{label}.stderr.txt"
                stdout_path.write_text(stdout_text, encoding="utf-8")
                stderr_path.write_text(stderr_text, encoding="utf-8")
            records.append(
                {
                    "label": label,
                    "command": command,
                    "returncode": returncode,
                    "duration_seconds": duration,
                    "stdout_path": str(stdout_path) if stdout_path is not None else None,
                    "stderr_path": str(stderr_path) if stderr_path is not None else None,
                    "stdout_tail": stdout_text.strip().splitlines()[-20:],
                    "stderr_tail": stderr_text.strip().splitlines()[-20:],
                }
            )
            if progress_callback is not None:
                progress_callback("finish", records[-1])
        return records
