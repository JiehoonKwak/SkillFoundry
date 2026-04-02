from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCRIPT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "openmm-forcefield-assignment-starter" / "scripts" / "run_openmm_forcefield_assignment.py"
INPUT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "openmm-forcefield-assignment-starter" / "examples" / "two_waters.pdb"
ASSET = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "openmm-forcefield-assignment-starter" / "assets" / "two_waters_forcefield_summary.json"


class OpenMmForceFieldAssignmentTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openmm_forcefield.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(INPUT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["residue_count"], 2)
            self.assertEqual(payload["atom_count"], 6)
            self.assertEqual(payload["particle_count"], 6)
            self.assertIn("NonbondedForce", payload["force_classes"])

    def test_committed_asset_matches_runtime_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openmm_forcefield.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(INPUT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertEqual(
                json.loads(out_path.read_text(encoding="utf-8")),
                json.loads(ASSET.read_text(encoding="utf-8")),
            )

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--input", "missing.pdb"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input PDB not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
