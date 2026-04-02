from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "reporting" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "reproducible-workflows"
    / "papermill-parameterized-notebook-starter"
    / "scripts"
    / "run_papermill_parameterized_notebook.py"
)
NOTEBOOK = (
    ROOT
    / "skills"
    / "reproducible-workflows"
    / "papermill-parameterized-notebook-starter"
    / "examples"
    / "toy_parameters.ipynb"
)


class PapermillParameterizedNotebookTests(unittest.TestCase):
    def test_executes_parameterized_notebook(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            output_notebook = tmp_root / "executed.ipynb"
            summary_out = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(NOTEBOOK),
                    "--output-notebook",
                    str(output_notebook),
                    "--summary-out",
                    str(summary_out),
                    "--x",
                    "5",
                    "--y",
                    "7",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertTrue(output_notebook.exists())
            self.assertEqual(payload["result"]["sum"], 12)
            self.assertEqual(payload["result"]["product"], 35)
            self.assertTrue(payload["injected_parameters_present"])
            self.assertGreaterEqual(payload["executed_code_cell_count"], 3)

    def test_rejects_missing_notebook(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--input",
                "missing.ipynb",
                "--output-notebook",
                "scratch/papermill/missing.ipynb",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=300,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Notebook not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
