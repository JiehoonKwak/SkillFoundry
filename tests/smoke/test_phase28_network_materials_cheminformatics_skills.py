from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
MATERIALS_PYTHON = ROOT / "slurm" / "envs" / "materials" / "bin" / "python"


class Phase28NetworkMaterialsCheminformaticsSmokeTests(unittest.TestCase):
    def test_rdkit_scaffold_analysis_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "scaffolds.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/drug-discovery-and-cheminformatics/rdkit-scaffold-analysis-starter/scripts/run_rdkit_scaffold_analysis.py",
                    "--input",
                    "skills/drug-discovery-and-cheminformatics/rdkit-scaffold-analysis-starter/examples/molecules.tsv",
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
            self.assertEqual(payload["summary"]["unique_murcko_scaffolds"], 2)
            self.assertEqual(payload["summary"]["largest_scaffold_group_size"], 3)

    def test_pymatgen_crystal_structure_parsing_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "structure.json"
            subprocess.run(
                [
                    str(MATERIALS_PYTHON),
                    "skills/materials-science-and-engineering/pymatgen-crystal-structure-parsing-starter/scripts/run_pymatgen_structure_summary.py",
                    "--input",
                    "skills/materials-science-and-engineering/pymatgen-crystal-structure-parsing-starter/examples/cscl.cif",
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
            self.assertEqual(payload["reduced_formula"], "CsCl")
            self.assertEqual(payload["symmetry"]["space_group_symbol"], "Pm-3m")

    def test_string_interaction_partners_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "string.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/string-interaction-partners-starter/scripts/run_string_interaction_partners.py",
                    "--identifier-file",
                    "skills/systems-biology/string-interaction-partners-starter/examples/tp53_query.txt",
                    "--species",
                    "9606",
                    "--limit",
                    "5",
                    "--required-score",
                    "700",
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
            self.assertEqual(payload["query_preferred_name"], "TP53")
            self.assertEqual(payload["partner_count"], 5)


if __name__ == "__main__":
    unittest.main()
