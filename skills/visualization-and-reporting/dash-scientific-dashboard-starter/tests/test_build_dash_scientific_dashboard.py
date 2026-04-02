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
    / "visualization-and-reporting"
    / "dash-scientific-dashboard-starter"
    / "scripts"
    / "build_dash_scientific_dashboard.py"
)
EXAMPLES = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "dash-scientific-dashboard-starter"
    / "examples"
)


class DashScientificDashboardTests(unittest.TestCase):
    def test_dashboard_outputs_are_generated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            html_out = tmp_root / "dashboard_preview.html"
            summary_out = tmp_root / "dashboard_summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLES / "toy_measurements.tsv"),
                    "--html-out",
                    str(html_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["measurement_count"], 5)
            self.assertEqual(payload["metric_options"], ["signal", "control"])
            self.assertEqual(payload["trace_count"], 2)
            self.assertEqual(payload["callback_count"], 1)
            self.assertIn("Dropdown", payload["layout_component_types"])
            self.assertTrue(Path(payload["shell_path"]).exists())

    def test_missing_input_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--input",
                "missing.tsv",
                "--html-out",
                "ignored.html",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Measurement table not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
