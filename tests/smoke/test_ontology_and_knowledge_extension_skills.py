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
        timeout=90,
    )
    return json.loads(completed.stdout)


class OntologyAndKnowledgeExtensionSkillSmokeTests(unittest.TestCase):
    def test_quickgo_term_search_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/systems-biology/quickgo-term-search/scripts/search_quickgo_terms.py",
            "--query",
            "apoptosis",
            "--limit",
            "1",
        )
        self.assertEqual(payload["query"], "apoptosis")
        self.assertEqual(payload["term_count"], 1)
        self.assertGreaterEqual(payload["number_of_hits"], 1)
        self.assertTrue(payload["terms"][0]["id"].startswith("GO:"))
        self.assertTrue(payload["terms"][0]["name"])


if __name__ == "__main__":
    unittest.main()
