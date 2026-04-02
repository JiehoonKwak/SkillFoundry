from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]


class BioSamplesSampleSearchTests(unittest.TestCase):
    def test_biosamples_search_returns_at_least_one_sample(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "biosamples.json"
            subprocess.run(
                [
                    "python3",
                    "skills/data-acquisition-and-dataset-handling/biosamples-sample-search/scripts/search_biosamples_samples.py",
                    "--text",
                    "breast cancer",
                    "--size",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["text"], "breast cancer")
            self.assertEqual(payload["sample_count"], 1)
            self.assertGreaterEqual(payload["total_elements"], 1)
            self.assertTrue(payload["samples"][0]["accession"].startswith("SAME"))
            self.assertTrue(payload["samples"][0]["sample_url"])


if __name__ == "__main__":
    unittest.main()
