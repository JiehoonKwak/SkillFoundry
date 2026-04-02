from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-computing-and-numerical-methods"
    / "scipy-ode-simulation-starter"
    / "scripts"
    / "run_scipy_ode_simulation.py"
)


class SciPyOdeSimulationTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_values(self) -> None:
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
            self.assertEqual(payload["model"], "lotka_volterra")
            self.assertAlmostEqual(payload["final_state"]["prey"], 0.537024, places=6)
            self.assertAlmostEqual(payload["final_state"]["predator"], 2.270535, places=6)
            self.assertAlmostEqual(payload["prey_peak_time"], 10.1, places=6)
            self.assertAlmostEqual(payload["predator_peak_time"], 0.8, places=6)

    def test_rejects_too_few_points(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--points", "8"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--points must be at least 20", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
