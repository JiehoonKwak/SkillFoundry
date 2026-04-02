from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "ase-geometry-optimization-starter" / "scripts" / "run_ase_geometry_optimization.py"
EXAMPLE = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "ase-geometry-optimization-starter" / "examples" / "toy_argon_dimer.json"
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"


class AseGeometryOptimizationTests(unittest.TestCase):
    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(CHEM_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertLess(payload["final_distance"], payload["initial_distance"])
            self.assertGreater(payload["energy_drop"], 0)
            self.assertAlmostEqual(payload["final_distance"], 3.81637, places=3)


if __name__ == "__main__":
    unittest.main()
