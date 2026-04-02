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
    / "gbif-dataset-search-starter"
    / "scripts"
    / "run_gbif_dataset_search.py"
)


class GbifDatasetSearchTests(unittest.TestCase):
    def test_runtime_summary_returns_dataset_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "datasets.json"
            subprocess.run(
                ["python3", str(SCRIPT), "--query", "puma", "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["query"], "puma")
            self.assertGreaterEqual(payload["result_count"], 1)
            self.assertIsNotNone(payload["first_dataset_key"])
            self.assertTrue(payload["first_title"])

    def test_rejects_invalid_limit(self) -> None:
        completed = subprocess.run(
            ["python3", str(SCRIPT), "--query", "puma", "--limit", "0"],
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
