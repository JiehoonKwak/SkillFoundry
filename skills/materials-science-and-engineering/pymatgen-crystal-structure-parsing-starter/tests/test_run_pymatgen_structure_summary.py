from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "materials" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "materials-science-and-engineering"
    / "pymatgen-crystal-structure-parsing-starter"
    / "scripts"
    / "run_pymatgen_structure_summary.py"
)
EXAMPLE = (
    ROOT
    / "skills"
    / "materials-science-and-engineering"
    / "pymatgen-crystal-structure-parsing-starter"
    / "examples"
    / "cscl.cif"
)
ASSET = (
    ROOT
    / "skills"
    / "materials-science-and-engineering"
    / "pymatgen-crystal-structure-parsing-starter"
    / "assets"
    / "cscl_structure_summary.json"
)


class PymatgenCrystalStructureParsingTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "structure_summary.json"
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
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["reduced_formula"], "CsCl")
            self.assertEqual(payload["site_count"], 2)
            self.assertEqual(payload["symmetry"]["space_group_symbol"], "Pm-3m")
            self.assertEqual(payload["symmetry"]["space_group_number"], 221)
            self.assertEqual(payload["species_counts"], {"Cl": 1, "Cs": 1})

    def test_committed_asset_matches_runtime_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "structure_summary.json"
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
                timeout=180,
            )
            self.assertEqual(
                json.loads(out_path.read_text(encoding="utf-8")),
                json.loads(ASSET.read_text(encoding="utf-8")),
            )

    def test_missing_input_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--input",
                "does-not-exist.cif",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("does not exist", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
