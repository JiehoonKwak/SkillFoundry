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
    / "mne-eeg-preprocessing-starter"
    / "scripts"
    / "run_mne_eeg_preprocessing.py"
)


class MneEegPreprocessingTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_filter_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "mne.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["channel_names"], ["Fz", "Cz"])
            self.assertEqual(payload["channel_count"], 2)
            self.assertAlmostEqual(payload["filtered_std_ratio"][0], 0.931382, places=6)
            self.assertAlmostEqual(payload["filtered_std_ratio"][1], 0.977195, places=6)


if __name__ == "__main__":
    unittest.main()
