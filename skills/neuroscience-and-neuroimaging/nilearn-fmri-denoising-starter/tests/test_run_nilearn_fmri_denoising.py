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
    / "nilearn-fmri-denoising-starter"
    / "scripts"
    / "run_nilearn_fmri_denoising.py"
)


class NilearnFmriDenoisingTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_confound_reduction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "nilearn.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["frame_count"], 60)
            self.assertEqual(payload["voxel_count"], 2)
            self.assertAlmostEqual(payload["pre_confound_abs_correlation_mean"], 0.384557, places=6)
            self.assertAlmostEqual(payload["post_confound_abs_correlation_mean"], 0.0, places=6)
            self.assertEqual(payload["cleaned_std_sample"], [1.0, 1.0])


if __name__ == "__main__":
    unittest.main()
