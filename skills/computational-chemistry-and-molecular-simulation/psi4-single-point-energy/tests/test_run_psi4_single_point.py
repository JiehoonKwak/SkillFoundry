from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "psi4-single-point-energy" / "scripts" / "run_psi4_single_point.py"
PSI4_PYTHON = ROOT / "slurm" / "envs" / "psi4" / "bin" / "python"


def load_script_module():
    spec = importlib.util.spec_from_file_location("psi4_single_point", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class Psi4SinglePointTests(unittest.TestCase):
    def test_build_molecule_metadata_returns_water_defaults(self) -> None:
        metadata = MODULE.build_molecule_metadata()
        self.assertEqual(metadata["name"], "water")
        self.assertEqual(metadata["formula"], "H2O")
        self.assertEqual(metadata["atom_count"], 3)
        self.assertEqual(metadata["geometry_angstrom"][0]["element"], "O")

    def test_runtime_single_point_summary(self) -> None:
        self.assertTrue(PSI4_PYTHON.exists(), f"Missing Psi4 interpreter: {PSI4_PYTHON}")
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "psi4-summary.json"
            subprocess.run(
                [str(PSI4_PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["mode"], "single_point")
            self.assertEqual(payload["method"], "hf")
            self.assertEqual(payload["basis"], "sto-3g")
            self.assertEqual(payload["molecule"]["atom_count"], 3)
            self.assertEqual(payload["molecule"]["formula"], "H2O")
            self.assertAlmostEqual(payload["energy_hartree"], -74.962991614813, places=9)
            self.assertEqual(payload["psi4_version"], "1.10")


if __name__ == "__main__":
    unittest.main()
