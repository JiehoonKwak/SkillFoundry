from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "causal" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "statistical-and-machine-learning-foundations-for-science"
    / "dowhy-average-treatment-effect-starter"
    / "scripts"
    / "run_dowhy_average_treatment_effect.py"
)


class DoWhyAverageTreatmentEffectTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_effect(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertAlmostEqual(payload["average_treatment_effect"], 1.855648, places=6)
            self.assertLess(abs(payload["placebo_effect"]), 0.1)
            self.assertGreater(payload["observed_outcome_gap"], payload["average_treatment_effect"])

    def test_rejects_too_small_sample_count(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--sample-count", "4"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("sample-count must be at least 10", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
