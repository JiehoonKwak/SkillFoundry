from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "statistical-and-machine-learning-foundations-for-science"
    / "optuna-bayesian-optimization-starter"
    / "scripts"
    / "run_optuna_bayesian_optimization.py"
)


class OptunaBayesianOptimizationTests(unittest.TestCase):
    def test_optimizes_toy_objective(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["trial_count"], 32)
            self.assertLess(payload["best_value"], 0.05)
            self.assertLess(abs(payload["best_params"]["x"] - 0.5), 0.2)
            self.assertLess(abs(payload["best_params"]["y"] + 0.25), 0.2)

    def test_rejects_too_few_trials(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--trials",
                "2",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("n_trials", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
