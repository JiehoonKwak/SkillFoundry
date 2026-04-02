from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "neuro" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "neuroscience-and-neuroimaging"
    / "nibabel-nifti-summary-starter"
    / "scripts"
    / "run_nibabel_nifti_summary.py"
)


class NibabelNiftiSummaryTests(unittest.TestCase):
    def test_generated_toy_nifti_summary_matches_expected_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            nifti_path = Path(tmp_dir) / "toy.nii.gz"
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--nifti-out", str(nifti_path), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertTrue(nifti_path.exists())
            self.assertEqual(payload["mode"], "generated_toy")
            self.assertEqual(payload["shape"], [4, 4, 3, 2])
            self.assertEqual(payload["zooms"], [2.0, 2.0, 2.0, 1.5])
            self.assertAlmostEqual(payload["mean_intensity"], 4.75, places=6)
            self.assertAlmostEqual(payload["affine_determinant"], 8.0, places=6)

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--input", "does-not-exist.nii.gz"],
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
