from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]


class EncodeExperimentSearchTests(unittest.TestCase):
    def test_encode_search_returns_at_least_one_result(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "encode.json"
            subprocess.run(
                [
                    "python3",
                    "skills/epigenomics-and-chromatin/encode-experiment-search/scripts/search_encode_experiments.py",
                    "--search-term",
                    "ATAC-seq",
                    "--limit",
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
            self.assertEqual(payload["search_term"], "ATAC-seq")
            self.assertEqual(payload["experiment_count"], 1)
            self.assertGreaterEqual(payload["total"], 1)
            self.assertTrue(payload["experiments"][0]["accession"])
            self.assertIn("ATAC-seq", payload["experiments"][0]["assay_title"])


if __name__ == "__main__":
    unittest.main()
