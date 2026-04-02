from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "deepchem" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "deepchem-circular-featurization"
    / "scripts"
    / "compute_circular_fingerprints.py"
)


class DeepChemCircularFeaturizationTests(unittest.TestCase):
    def test_default_example_writes_expected_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fingerprints.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["featurizer"], "CircularFingerprint")
            self.assertEqual(payload["deepchem_version"], "2.8.0")
            self.assertEqual(payload["size"], 32)
            self.assertEqual(payload["radius"], 2)
            self.assertEqual(payload["molecule_count"], 2)
            self.assertEqual(payload["molecules"][0]["canonical_smiles"], "CC(=O)Oc1ccccc1C(=O)O")
            self.assertEqual(len(payload["molecules"][0]["bit_vector"]), 32)
            self.assertEqual(payload["molecules"][0]["on_bit_count"], len(payload["molecules"][0]["on_bits"]))

    def test_single_smiles_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            first = Path(tmp_dir) / "first.json"
            second = Path(tmp_dir) / "second.json"
            command = [str(PYTHON), str(SCRIPT), "--smiles", "CCO", "--size", "16", "--out"]
            subprocess.run(command + [str(first)], cwd=ROOT, check=True, capture_output=True, text=True, timeout=60)
            subprocess.run(command + [str(second)], cwd=ROOT, check=True, capture_output=True, text=True, timeout=60)
            self.assertEqual(
                json.loads(first.read_text(encoding="utf-8")),
                json.loads(second.read_text(encoding="utf-8")),
            )

    def test_invalid_smiles_fails_cleanly(self) -> None:
        result = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--smiles", "not-a-smiles"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=60,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid SMILES", result.stderr or result.stdout)


if __name__ == "__main__":
    unittest.main()
