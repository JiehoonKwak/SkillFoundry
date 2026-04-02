from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "astronomy" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "physics-and-astronomy"
    / "astropy-fits-image-summary-starter"
    / "scripts"
    / "run_astropy_fits_image_summary.py"
)


class AstropyFitsImageSummaryTests(unittest.TestCase):
    def test_generated_toy_fits_summary_matches_expected_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            fits_path = Path(tmp_dir) / "toy.fits"
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--fits-out", str(fits_path), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertTrue(fits_path.exists())
            self.assertEqual(payload["shape"], [5, 5])
            self.assertAlmostEqual(payload["mean_signal"], 12.0, places=6)
            self.assertAlmostEqual(payload["std_signal"], 7.211102, places=5)
            self.assertEqual(payload["filter"], "r")
            self.assertAlmostEqual(payload["exptime_s"], 30.0, places=6)

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--input", "missing.fits"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input file not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
