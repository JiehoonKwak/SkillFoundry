from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "metagenomics" / "bin" / "python"
SCRIPT = ROOT / "skills" / "genomics" / "sourmash-signature-compare-starter" / "scripts" / "run_sourmash_signature_compare.py"


class SourmashSignatureCompareTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_similarity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["shared_hash_count"], 12)
            self.assertAlmostEqual(payload["jaccard_similarity"], 0.631579, places=6)
            self.assertAlmostEqual(payload["query_containment_in_reference"], 0.8, places=6)

    def test_rejects_invalid_scaled(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--scaled", "0"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("scaled must be a positive integer", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
