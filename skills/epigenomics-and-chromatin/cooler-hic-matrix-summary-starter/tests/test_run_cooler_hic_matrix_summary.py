from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
SCRIPT = ROOT / "skills" / "epigenomics-and-chromatin" / "cooler-hic-matrix-summary-starter" / "scripts" / "run_cooler_hic_matrix_summary.py"


class CoolerHicMatrixSummaryTests(unittest.TestCase):
    def test_generated_toy_summary_matches_expected_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "cooler_summary.json"
            cool_path = Path(tmp_dir) / "toy.cool"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--cooler-out",
                    str(cool_path),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["bin_count"], 3)
            self.assertEqual(payload["pixel_count"], 5)
            self.assertEqual(payload["matrix_shape"], [3, 3])
            self.assertEqual(payload["total_contact_count"], 38)
            self.assertEqual(payload["first_pixels"][0]["count"], 10)
            self.assertTrue(cool_path.exists())


if __name__ == "__main__":
    unittest.main()
