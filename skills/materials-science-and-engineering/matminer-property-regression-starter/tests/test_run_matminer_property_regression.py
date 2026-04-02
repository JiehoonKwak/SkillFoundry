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
    / "matminer-property-regression-starter"
    / "scripts"
    / "run_matminer_property_regression.py"
)


class MatminerPropertyRegressionTests(unittest.TestCase):
    def test_runtime_summary_contains_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "materials_summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
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
            self.assertGreater(payload["feature_count"], 100)
            self.assertEqual(payload["test_count"], 2)
            self.assertEqual(len(payload["holdout_predictions"]), 2)


if __name__ == "__main__":
    unittest.main()
