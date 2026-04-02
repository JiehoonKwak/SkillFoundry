from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "frictionless-tabular-validation-starter"
    / "scripts"
    / "run_frictionless_tabular_validation.py"
)
EXAMPLES = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "frictionless-tabular-validation-starter"
    / "examples"
)


class FrictionlessTabularValidationTests(unittest.TestCase):
    def test_valid_table_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLES / "toy_people_valid.csv"),
                    "--schema",
                    str(EXAMPLES / "toy_people_schema.json"),
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
            self.assertTrue(payload["valid"])
            self.assertEqual(payload["row_count"], 3)
            self.assertEqual(payload["error_count"], 0)
            self.assertEqual(payload["field_names"], ["sample_id", "condition", "count"])

    def test_invalid_table_reports_type_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLES / "toy_people_invalid.csv"),
                    "--schema",
                    str(EXAMPLES / "toy_people_schema.json"),
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
            self.assertFalse(payload["valid"])
            self.assertEqual(payload["error_count"], 1)
            self.assertEqual(payload["errors"][0]["type"], "type-error")

    def test_missing_input_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--input",
                "missing.csv",
                "--schema",
                str(EXAMPLES / "toy_people_schema.json"),
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input table not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
