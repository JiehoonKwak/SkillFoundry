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
    / "drug-discovery-and-cheminformatics"
    / "rdkit-molecule-standardization"
    / "scripts"
    / "standardize_rdkit_molecule.py"
)
EXAMPLE = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "rdkit-molecule-standardization"
    / "examples"
    / "sodium_acetate.smiles"
)
ASSET = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "rdkit-molecule-standardization"
    / "assets"
    / "sodium_acetate_standardized.json"
)


class RdkitMoleculeStandardizationTests(unittest.TestCase):
    def test_runtime_standardization_strips_salt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--smiles-file",
                    str(EXAMPLE),
                    "--name",
                    "sodium-acetate",
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
            self.assertEqual(payload["canonical_tautomer_smiles"], "CC(=O)O")
            self.assertTrue(payload["removed_fragments"])
            self.assertEqual(payload["charge_after"], 0)

    def test_committed_asset_matches_runtime_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--smiles-file",
                    str(EXAMPLE),
                    "--name",
                    "sodium-acetate",
                    "--out",
                    str(out_path),
                ],
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

    def test_invalid_smiles_returns_error(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--smiles", "not-a-smiles"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Invalid SMILES", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
