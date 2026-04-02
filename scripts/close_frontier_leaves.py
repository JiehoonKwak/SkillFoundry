#!/usr/bin/env python3
"""Sync missing skills and generate starter skills for all frontier leaves."""

from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path

from build_site import TAXONOMY_DOMAIN_MAP, build_tree, load_json, load_jsonl


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "registry"
SKILLS = ROOT / "skills"
TODAY = date.today().isoformat()
SAFE_SLUG_RE = re.compile(r"[^a-z0-9]+")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def safe_slug(value: str) -> str:
    return SAFE_SLUG_RE.sub("-", value.lower()).strip("-")


def load_registry_resources() -> list[dict]:
    return load_jsonl(REGISTRY / "resources_dedup.jsonl")


def load_registry_skills() -> list[dict]:
    return load_jsonl(REGISTRY / "skills.jsonl")


def registry_entry_from_metadata(metadata: dict, skill_dir: Path) -> dict:
    source_resource_ids = metadata.get("source_resource_ids") or metadata.get("source_resources") or []
    return {
        "skill_id": metadata.get("skill_id", metadata["slug"]),
        "name": metadata["name"],
        "slug": metadata["slug"],
        "domain": metadata["domain"],
        "topic_path": metadata["topic_path"],
        "status": metadata["status"],
        "skill_type": metadata.get("skill_type", "starter"),
        "description": metadata.get("description", f"Starter skill for {metadata['name'].lower()}."),
        "path": str(skill_dir.relative_to(ROOT)),
        "source_resource_ids": source_resource_ids,
        "tags": metadata.get("tags", []),
        "compute_requirements": metadata.get("compute_requirements", "cpu+local"),
        "language": metadata.get("language", ["Python"]),
        "test_commands": metadata.get("test_commands", []),
        "last_updated": metadata.get("updated_at", metadata.get("last_updated", TODAY)),
        "last_verified": metadata.get("last_verified_at", metadata.get("last_verified", TODAY)),
        "related_skills": metadata.get("related_skills", []),
    }


def discover_missing_skill_entries(existing_slugs: set[str]) -> list[dict]:
    discovered: list[dict] = []
    for metadata_path in sorted(SKILLS.glob("*/*/metadata.yaml")):
        metadata = read_json(metadata_path)
        slug = metadata["slug"]
        if slug in existing_slugs:
            continue
        discovered.append(registry_entry_from_metadata(metadata, metadata_path.parent))
        existing_slugs.add(slug)
    return discovered


def frontier_leaf_nodes(skills: list[dict], resources: list[dict], taxonomy: dict) -> list[dict]:
    tree = build_tree(skills, resources, taxonomy)
    leaves: list[dict] = []
    for domain in tree["children"]:
        domain_slug = TAXONOMY_DOMAIN_MAP[domain["taxonomy_key"]][0]
        for leaf in domain["children"]:
            if leaf["coverage_status"] != "frontier":
                continue
            leaves.append(
                {
                    "domain_name": domain["name"],
                    "domain_slug": domain_slug,
                    "taxonomy_key": domain["taxonomy_key"],
                    "leaf_name": leaf["name"],
                    "leaf_slug": leaf["topic_slug"],
                    "resources": leaf["resources"],
                    "resource_ids": [resource["resource_id"] for resource in leaf["resources"]],
                }
            )
    return leaves


def skill_metadata_payload(leaf: dict, skill_slug: str) -> dict:
    skill_name = f"{leaf['leaf_name']} Starter"
    example_out = f"scratch/frontier/{skill_slug}.json"
    return {
        "name": skill_name,
        "slug": skill_slug,
        "domain": leaf["domain_slug"],
        "topic_path": [leaf["domain_slug"], leaf["leaf_slug"], "starter"],
        "status": "implemented",
        "skill_type": "frontier-starter",
        "description": f"Generate an actionable starter plan for the {leaf['leaf_name']} frontier leaf from curated local resources.",
        "authorship": "Codex",
        "created_at": TODAY,
        "updated_at": TODAY,
        "source_resources": leaf["resource_ids"],
        "source_resource_ids": leaf["resource_ids"],
        "package_versions": {},
        "compute_requirements": "cpu+local",
        "language": ["Python", "JSON", "Markdown"],
        "tags": ["frontier-closure", leaf["domain_slug"], leaf["leaf_slug"]],
        "aliases": [
            leaf["leaf_name"].lower(),
            leaf["leaf_slug"].replace("-", " "),
            skill_slug.replace("-", " "),
        ],
        "test_commands": [
            f"python3 skills/{leaf['domain_slug']}/{skill_slug}/scripts/run_frontier_starter.py --out {example_out}"
        ],
        "last_verified_at": TODAY,
        "freshness_risk": "low",
        "related_skills": [],
    }


def skill_registry_entry(leaf: dict, skill_slug: str) -> dict:
    return {
        "skill_id": skill_slug,
        "name": f"{leaf['leaf_name']} Starter",
        "slug": skill_slug,
        "domain": leaf["domain_slug"],
        "topic_path": [leaf["domain_slug"], leaf["leaf_slug"], "starter"],
        "status": "implemented",
        "skill_type": "frontier-starter",
        "description": f"Starter plan generator for the {leaf['leaf_name']} frontier leaf.",
        "path": f"skills/{leaf['domain_slug']}/{skill_slug}",
        "source_resource_ids": leaf["resource_ids"],
        "tags": ["frontier-closure", leaf["domain_slug"], leaf["leaf_slug"]],
        "compute_requirements": "cpu+local",
        "language": ["Python", "JSON", "Markdown"],
        "test_commands": [
            f"python3 skills/{leaf['domain_slug']}/{skill_slug}/scripts/run_frontier_starter.py --out scratch/frontier/{skill_slug}.json"
        ],
        "last_updated": TODAY,
        "last_verified": TODAY,
        "related_skills": [],
    }


def skill_files(leaf: dict, skill_slug: str) -> dict[str, str]:
    skill_name = f"{leaf['leaf_name']} Starter"
    refs_lines = [
        "# References",
        "",
        f"This starter closes the `{leaf['leaf_slug']}` frontier leaf in `{leaf['domain_name']}`.",
        "",
    ]
    for resource in leaf["resources"]:
        refs_lines.append(
            f"- `{resource['resource_id']}`: [{resource['canonical_name']}]({resource['url']})"
        )
    refs_lines.append("")
    refs_lines.append("Use these curated anchors to promote this starter into a runtime-verified skill.")
    example_context = {
        "skill_slug": skill_slug,
        "skill_name": skill_name,
        "domain_name": leaf["domain_name"],
        "domain_slug": leaf["domain_slug"],
        "taxonomy_key": leaf["taxonomy_key"],
        "leaf_name": leaf["leaf_name"],
        "leaf_slug": leaf["leaf_slug"],
        "source_resource_ids": leaf["resource_ids"],
        "starter_objectives": [
            f"Review the primary materials for {leaf['leaf_name']}.",
            "Define the smallest reproducible input/output contract.",
            "Capture a smoke command or toy example.",
            "Promote the starter to sandbox verification once runtime details are stable.",
        ],
    }
    script = f"""#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=SKILL_ROOT / "assets" / "starter_summary.json")
    args = parser.parse_args()
    metadata = json.loads((SKILL_ROOT / "metadata.yaml").read_text(encoding="utf-8"))
    context = json.loads((SKILL_ROOT / "examples" / "resource_context.json").read_text(encoding="utf-8"))
    summary = {{
        "skill_slug": metadata["slug"],
        "skill_name": metadata["name"],
        "status": metadata["status"],
        "leaf_slug": context["leaf_slug"],
        "leaf_name": context["leaf_name"],
        "domain_slug": context["domain_slug"],
        "source_resource_ids": context["source_resource_ids"],
        "starter_steps": context["starter_objectives"],
        "promotion_checklist": [
            "Add a runnable example or toy dataset.",
            "Add a repository-level smoke or integration test.",
            "Promote status to sandbox_verified after checks pass."
        ]
    }}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
"""
    test_script = f"""from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "{leaf['domain_slug']}" / "{skill_slug}"


class FrontierStarterTests(unittest.TestCase):
    def test_starter_summary_contains_leaf_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "summary.json"
            subprocess.run(
                ["python3", str(SKILL / "scripts" / "run_frontier_starter.py"), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["skill_slug"], "{skill_slug}")
            self.assertEqual(payload["leaf_slug"], "{leaf['leaf_slug']}")
            self.assertEqual(payload["source_resource_ids"], {json.dumps(leaf["resource_ids"])})


if __name__ == "__main__":
    unittest.main()
"""
    return {
        "SKILL.md": "\n".join(
            [
                f"# {skill_name}",
                "",
                f"Use this starter when a task lands in the `{leaf['leaf_name']}` frontier leaf and the repository has curated resources but no dedicated runtime implementation yet.",
                "",
                "## What this starter does",
                "",
                "- Summarizes the local resource anchors for the leaf.",
                "- Emits a machine-readable starter plan with promotion steps.",
                "- Gives the agent a stable local entry point before a full runtime skill exists.",
                "",
                "## How to use it",
                "",
                f"Run `python3 skills/{leaf['domain_slug']}/{skill_slug}/scripts/run_frontier_starter.py --out scratch/frontier/{skill_slug}.json`.",
                "",
                "Then inspect `refs.md` and `examples/resource_context.json` to promote the starter into a concrete executable workflow.",
                "",
            ]
        )
        + "\n",
        "metadata.yaml": json.dumps(skill_metadata_payload(leaf, skill_slug), indent=2) + "\n",
        "refs.md": "\n".join(refs_lines) + "\n",
        "examples/README.md": (
            "# Example Context\n\n"
            "The generated `resource_context.json` file captures the exact taxonomy leaf and source resources that back this starter.\n"
        ),
        "examples/resource_context.json": json.dumps(example_context, indent=2) + "\n",
        "assets/README.md": "# Assets\n\nGenerated starter summaries can be written here or into `scratch/frontier/` during smoke runs.\n",
        "scripts/run_frontier_starter.py": script,
        f"tests/test_{skill_slug.replace('-', '_')}.py": test_script,
    }


def write_skill_tree(skill_dir: Path, files: dict[str, str]) -> None:
    for relative_path, content in files.items():
        path = skill_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")


def append_aliases(existing: dict[str, list[str]], entry: dict) -> None:
    key = entry["slug"].replace("-", "_")
    phrases = [
        entry["name"].lower(),
        entry["slug"].replace("-", " "),
        entry["topic_path"][1].replace("-", " "),
    ]
    current = existing.setdefault(key, [])
    for phrase in phrases:
        if phrase not in current:
            current.append(phrase)


def main() -> int:
    resources = load_registry_resources()
    skills = load_registry_skills()
    taxonomy = load_json(REGISTRY / "taxonomy.yaml")
    aliases = read_json(REGISTRY / "aliases.json")

    existing_slugs = {skill["slug"] for skill in skills}
    synced_entries = discover_missing_skill_entries(existing_slugs)
    for entry in synced_entries:
        skills.append(entry)
        append_aliases(aliases, entry)

    leaves = frontier_leaf_nodes(skills, resources, taxonomy)
    generated_entries: list[dict] = []
    for leaf in leaves:
        skill_slug = f"{safe_slug(leaf['leaf_slug'])}-starter"
        if skill_slug in existing_slugs:
            continue
        skill_dir = SKILLS / leaf["domain_slug"] / skill_slug
        write_skill_tree(skill_dir, skill_files(leaf, skill_slug))
        entry = skill_registry_entry(leaf, skill_slug)
        skills.append(entry)
        generated_entries.append(entry)
        existing_slugs.add(skill_slug)
        append_aliases(aliases, entry)

    (REGISTRY / "skills.jsonl").write_text(
        "".join(json.dumps(skill) + "\n" for skill in skills),
        encoding="utf-8",
    )
    write_json(REGISTRY / "aliases.json", aliases)

    updated_tree = build_tree(skills, resources, taxonomy)
    print(
        json.dumps(
            {
                "synced_existing_skills": len(synced_entries),
                "generated_skills": len(generated_entries),
                "skill_count": updated_tree["skill_count"],
                "covered_leaf_count": updated_tree["covered_leaf_count"],
                "frontier_leaf_count": updated_tree["frontier_leaf_count"],
                "todo_leaf_count": updated_tree["todo_leaf_count"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
