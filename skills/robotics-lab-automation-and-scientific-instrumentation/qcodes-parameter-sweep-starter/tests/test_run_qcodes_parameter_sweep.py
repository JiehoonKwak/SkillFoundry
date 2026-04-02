from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "instrumentation" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "robotics-lab-automation-and-scientific-instrumentation"
    / "qcodes-parameter-sweep-starter"
    / "scripts"
    / "run_qcodes_parameter_sweep.py"
)
SETPOINTS = (
    ROOT
    / "skills"
    / "robotics-lab-automation-and-scientific-instrumentation"
    / "qcodes-parameter-sweep-starter"
    / "examples"
    / "toy_setpoints.tsv"
)


class QCoDeSParameterSweepTests(unittest.TestCase):
    def test_runs_parameter_sweep_and_records_dataset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            db_out = tmp_root / "toy_sweep.db"
            out_path = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--setpoints",
                    str(SETPOINTS),
                    "--db-out",
                    str(db_out),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["run_id"], 1)
            self.assertEqual(payload["point_count"], 4)
            self.assertEqual(payload["setpoint_min"], 0.0)
            self.assertEqual(payload["setpoint_max"], 1.5)
            self.assertAlmostEqual(payload["response_mean"], 1.15, places=6)
            self.assertTrue(db_out.exists())

    def test_missing_setpoints_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--setpoints",
                "missing.tsv",
                "--db-out",
                "scratch/qcodes/missing.db",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Setpoints TSV not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
