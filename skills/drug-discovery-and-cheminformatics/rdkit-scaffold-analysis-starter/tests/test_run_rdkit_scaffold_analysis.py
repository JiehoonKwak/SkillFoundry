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
    / "rdkit-scaffold-analysis-starter"
    / "scripts"
    / "run_rdkit_scaffold_analysis.py"
)
EXAMPLE = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "rdkit-scaffold-analysis-starter"
    / "examples"
    / "molecules.tsv"
)
ASSET = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "rdkit-scaffold-analysis-starter"
    / "assets"
    / "example_scaffold_summary.json"
)


class RdkitScaffoldAnalysisTests(unittest.TestCase):
    def test_example_summary_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "scaffolds.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLE),
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
            self.assertEqual(payload["summary"]["molecule_count"], 4)
            self.assertEqual(payload["summary"]["unique_murcko_scaffolds"], 2)
            self.assertEqual(payload["summary"]["unique_generic_scaffolds"], 2)
            self.assertEqual(payload["summary"]["largest_scaffold_group_size"], 3)
            self.assertEqual(payload["scaffold_groups"][0]["murcko_scaffold"], "c1ccccc1")
            self.assertEqual(payload["scaffold_groups"][0]["count"], 3)

    def test_committed_asset_matches_runtime_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "scaffolds.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLE),
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
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_input = Path(tmp_dir) / "bad.tsv"
            bad_input.write_text(
                "name\tsmiles\nbad-entry\tnot_a_smiles\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(bad_input),
                ],
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
