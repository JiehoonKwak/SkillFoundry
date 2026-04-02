from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "quarto-notebook-report-starter"
    / "scripts"
    / "render_quarto_notebook_report.py"
)
NOTEBOOK = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "quarto-notebook-report-starter"
    / "examples"
    / "toy_report.ipynb"
)


class QuartoNotebookReportTests(unittest.TestCase):
    def test_renders_html_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            html_out = Path(tmp_dir) / "report.html"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    str(NOTEBOOK),
                    "--html-out",
                    str(html_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertTrue(html_out.exists())
            self.assertTrue(payload["title_present"])
            self.assertTrue(payload["executed_output_present"])
            self.assertGreater(payload["html_size_bytes"], 1000)

    def test_rejects_missing_notebook(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                "missing-notebook.ipynb",
                "--html-out",
                "scratch/quarto/missing.html",
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
