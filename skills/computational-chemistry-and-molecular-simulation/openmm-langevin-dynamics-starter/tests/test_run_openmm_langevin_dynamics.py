from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "computational-chemistry-and-molecular-simulation"
    / "openmm-langevin-dynamics-starter"
    / "scripts"
    / "run_openmm_langevin_dynamics.py"
)


class OpenmmLangevinDynamicsTests(unittest.TestCase):
    def test_generates_deterministic_summary(self) -> None:
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
            self.assertEqual(payload["frame_count"], 6)
            self.assertEqual(payload["frames"][0]["step"], 0)
            self.assertEqual(payload["frames"][-1]["step"], 50)
            self.assertGreater(payload["final_bond_length_nm"], payload["initial_bond_length_nm"])
            self.assertLess(payload["frames"][0]["potential_energy_kj_mol"], payload["frames"][-1]["potential_energy_kj_mol"])
            self.assertLess(payload["min_potential_energy_kj_mol"], payload["max_potential_energy_kj_mol"])

    def test_rejects_bad_report_interval(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--steps",
                "10",
                "--report-interval",
                "20",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("report_interval", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
