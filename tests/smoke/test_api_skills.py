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
        timeout=60,
    )
    return json.loads(completed.stdout)


class ApiSkillSmokeTests(unittest.TestCase):
    def test_openalex_search_returns_results(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/openalex-literature-search/scripts/search_openalex.py",
            "--query",
            "single-cell RNA-seq",
            "--per-page",
            "1",
        )
        self.assertIn("meta", payload)
        self.assertIn("results", payload)
        self.assertGreaterEqual(len(payload["results"]), 1)

    def test_europepmc_search_returns_results(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/europepmc-method-triage/scripts/search_europepmc.py",
            "--query",
            "single-cell RNA-seq",
            "--page-size",
            "1",
        )
        self.assertIn("hitCount", payload)
        self.assertIn("resultList", payload)
        self.assertGreaterEqual(len(payload["resultList"]["result"]), 1)


if __name__ == "__main__":
    unittest.main()
