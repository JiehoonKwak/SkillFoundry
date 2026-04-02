from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCRIPT = ROOT / "skills" / "drug-discovery-and-cheminformatics" / "rdkit-molecular-descriptors" / "scripts" / "compute_rdkit_descriptors.py"


class RdkitDescriptorTests(unittest.TestCase):
    def test_aspirin_descriptor_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "aspirin.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--smiles",
                    "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "--name",
                    "aspirin",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["name"], "aspirin")
            self.assertEqual(payload["formula"], "C9H8O4")
            self.assertGreater(payload["molecular_weight"], 180.0)
            self.assertEqual(payload["hbd"], 1)


if __name__ == "__main__":
    unittest.main()

