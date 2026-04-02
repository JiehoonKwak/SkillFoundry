from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "ecology-evolution-and-biodiversity"
    / "gbif-species-occurrence-search-starter"
    / "scripts"
    / "run_gbif_species_occurrence_search.py"
)


class GbifSpeciesOccurrenceSearchTests(unittest.TestCase):
    def test_runtime_summary_returns_species_match_and_occurrences(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "gbif.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--scientific-name",
                    "Puma concolor",
                    "--country",
                    "US",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["query_scientific_name"], "Puma concolor")
            self.assertEqual(payload["country"], "US")
            self.assertGreaterEqual(payload["occurrence_count"], 1)
            self.assertTrue(payload["matched_scientific_name"].startswith("Puma concolor"))

    def test_rejects_invalid_limit(self) -> None:
        completed = subprocess.run(
            ["python3", str(SCRIPT), "--scientific-name", "Puma concolor", "--limit", "0"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--limit must be between 1 and 20", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
