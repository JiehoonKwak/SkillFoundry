from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCRIPT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "openmm-system-minimization" / "scripts" / "run_openmm_minimization.py"


class OpenMmMinimizationTests(unittest.TestCase):
    def test_energy_decreases_after_minimization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openmm.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["particle_count"], 3)
            self.assertLess(payload["final_potential_energy_kj_mol"], payload["initial_potential_energy_kj_mol"])
            self.assertEqual(len(payload["final_positions_nm"]), 3)


if __name__ == "__main__":
    unittest.main()
