from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run_json_command(*args: str) -> dict:
    completed = subprocess.run(
        ["python3", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
    )
    return json.loads(completed.stdout)


class CompositeKnowledgeSkillSmokeTests(unittest.TestCase):
    def test_multi_source_literature_brief_returns_all_sources(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/multi-source-literature-brief/scripts/build_literature_brief.py",
            "--query",
            "single-cell RNA-seq",
            "--limit",
            "1",
        )
        self.assertEqual(payload["limit_per_source"], 1)
        self.assertEqual(sorted(payload["sources"].keys()), ["crossref", "europepmc", "openalex", "pubmed"])
        self.assertTrue(all(payload["sources"][key] for key in payload["sources"]))


if __name__ == "__main__":
    unittest.main()
