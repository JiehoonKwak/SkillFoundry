from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCRIPT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "rdkit-conformer-generation-starter" / "scripts" / "run_rdkit_conformer_generation.py"
INPUT = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "rdkit-conformer-generation-starter" / "examples" / "molecules.tsv"
ASSET = ROOT / "skills" / "computational-chemistry-and-molecular-simulation" / "rdkit-conformer-generation-starter" / "assets" / "example_conformer_summary.json"


class RdkitConformerGenerationTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "conformers.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(INPUT), "--num-confs", "4", "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["molecule_count"], 2)
            self.assertEqual(payload["num_confs_requested"], 4)
            self.assertEqual(payload["molecules"][0]["conformer_count"], 4)
            self.assertEqual(payload["molecules"][1]["conformer_count"], 4)

    def test_committed_asset_matches_runtime_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "conformers.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(INPUT), "--num-confs", "4", "--out", str(out_path)],
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
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_input = Path(tmp_dir) / "bad.tsv"
            bad_input.write_text("name\tsmiles\nbad\tnot_a_smiles\n", encoding="utf-8")
            completed = subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(bad_input)],
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
