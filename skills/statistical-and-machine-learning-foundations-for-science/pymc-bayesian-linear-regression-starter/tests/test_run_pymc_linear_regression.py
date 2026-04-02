from __future__ import annotations

import importlib.util
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
    / "pymc-bayesian-linear-regression-starter"
    / "scripts"
    / "run_pymc_linear_regression.py"
)
INPUT = (
    ROOT
    / "skills"
    / "statistical-and-machine-learning-foundations-for-science"
    / "pymc-bayesian-linear-regression-starter"
    / "examples"
    / "toy_observations.tsv"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("pymc_linear_regression", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class PyMCLinearRegressionTests(unittest.TestCase):
    def test_module_exposes_expected_entry_points(self) -> None:
        self.assertTrue(hasattr(MODULE, "load_observations"))
        self.assertTrue(hasattr(MODULE, "fit_linear_regression"))

    def test_runtime_linear_regression_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "pymc-summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--draws",
                    "120",
                    "--tune",
                    "120",
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
            self.assertGreater(payload["posterior_means"]["slope"], 0.7)
            self.assertLess(payload["posterior_means"]["slope"], 1.0)
            self.assertGreater(payload["posterior_means"]["intercept"], 1.0)
            self.assertLess(payload["posterior_means"]["intercept"], 1.6)
            self.assertGreater(payload["hdi_90"]["slope"][0], 0.75)
            self.assertLess(payload["hdi_90"]["slope"][1], 0.9)
            self.assertGreater(payload["posterior_predictive_mean"][-1]["mean_y"], 5.0)


if __name__ == "__main__":
    unittest.main()
