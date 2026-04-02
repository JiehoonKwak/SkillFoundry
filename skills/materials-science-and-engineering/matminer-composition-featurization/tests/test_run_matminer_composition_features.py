from __future__ import annotations

import importlib.util
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
    / "matminer-composition-featurization"
    / "scripts"
    / "run_matminer_composition_features.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("matminer_composition_features", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class MatminerCompositionFeatureTests(unittest.TestCase):
    def test_normalize_formulas_uses_defaults(self) -> None:
        self.assertEqual(MODULE.normalize_formulas(None), ["Fe2O3", "LiFePO4", "SiO2"])

    def test_runtime_composition_features(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "matminer.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--formula", "Fe2O3", "--formula", "LiFePO4", "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["formula_count"], 2)
            self.assertEqual(payload["summaries"][0]["reduced_formula"], "Fe2O3")
            self.assertEqual(payload["summaries"][0]["stoichiometry_features"]["0-norm"], 2.0)
            self.assertEqual(payload["summaries"][1]["reduced_formula"], "LiFePO4")


if __name__ == "__main__":
    unittest.main()
