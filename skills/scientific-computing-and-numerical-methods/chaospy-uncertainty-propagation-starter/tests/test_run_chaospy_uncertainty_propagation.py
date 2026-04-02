from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "numerics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-computing-and-numerical-methods"
    / "chaospy-uncertainty-propagation-starter"
    / "scripts"
    / "run_chaospy_uncertainty_propagation.py"
)
EXAMPLES = (
    ROOT
    / "skills"
    / "scientific-computing-and-numerical-methods"
    / "chaospy-uncertainty-propagation-starter"
    / "examples"
)


class ChaospyUncertaintyPropagationTests(unittest.TestCase):
    def test_summary_matches_expected_distribution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--config",
                    str(EXAMPLES / "toy_parameters.json"),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["node_count"], 16)
            self.assertAlmostEqual(payload["mean"], 1.51, places=6)
            self.assertAlmostEqual(payload["std"], 0.202567, places=6)
            self.assertLess(payload["percentiles"]["p05"], payload["percentiles"]["p50"])
            self.assertLess(payload["percentiles"]["p50"], payload["percentiles"]["p95"])

    def test_missing_config_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--config",
                "missing.json",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Config not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
