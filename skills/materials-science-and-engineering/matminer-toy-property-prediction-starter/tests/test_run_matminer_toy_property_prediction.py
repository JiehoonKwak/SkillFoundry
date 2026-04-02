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
    / "matminer-toy-property-prediction-starter"
    / "scripts"
    / "run_matminer_toy_property_prediction.py"
)
INPUT = (
    ROOT
    / "skills"
    / "materials-science-and-engineering"
    / "matminer-toy-property-prediction-starter"
    / "examples"
    / "toy_materials.tsv"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("run_matminer_toy_property_prediction", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class MatminerToyPropertyPredictionTests(unittest.TestCase):
    def test_rank_feature_importances_orders_descending(self) -> None:
        ranked = MODULE.rank_feature_importances(["b", "a", "c"], [0.2, 0.5, 0.0], limit=2)
        self.assertEqual(ranked[0]["feature"], "a")
        self.assertEqual(ranked[1]["feature"], "b")

    def test_runtime_prediction_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "materials_summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
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
            self.assertEqual(payload["sample_count"], 6)
            self.assertEqual(len(payload["predictions"]), 6)
            self.assertLessEqual(payload["training_mae"], 0.001)
            self.assertGreaterEqual(len(payload["top_feature_importances"]), 1)


if __name__ == "__main__":
    unittest.main()
