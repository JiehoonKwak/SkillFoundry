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


class PathwayApiSkillSmokeTests(unittest.TestCase):
    def test_reactome_event_summary_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/systems-biology/reactome-event-summary/scripts/fetch_reactome_event_summary.py",
            "--stable-id",
            "R-HSA-141409",
        )
        self.assertEqual(payload["stable_id"], "R-HSA-141409")
        self.assertEqual(payload["species"], "Homo sapiens")
        self.assertTrue(payload["display_name"])
        self.assertTrue(payload["source_url"].endswith("/R-HSA-141409"))

    def test_reactome_identifiers_enrichment_returns_pathways(self) -> None:
        payload = run_json_command(
            "skills/systems-biology/reactome-identifiers-enrichment/scripts/analyze_reactome_identifiers.py",
            "--identifiers",
            "BRCA1,TP53",
            "--page-size",
            "3",
        )
        self.assertEqual(payload["identifiers"], ["BRCA1", "TP53"])
        self.assertGreaterEqual(len(payload["pathways"]), 1)
        self.assertTrue(payload["pathways"][0]["stable_id"].startswith("R-HSA-"))
        self.assertEqual(payload["pathways"][0]["species"], "Homo sapiens")
        self.assertIsNotNone(payload["pathways"][0]["entities_fdr"])


if __name__ == "__main__":
    unittest.main()
