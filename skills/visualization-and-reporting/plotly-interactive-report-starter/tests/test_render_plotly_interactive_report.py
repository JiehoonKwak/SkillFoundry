from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "plotly-interactive-report-starter"
    / "scripts"
    / "render_plotly_interactive_report.py"
)
INPUT = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "plotly-interactive-report-starter"
    / "examples"
    / "toy_measurements.tsv"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("render_plotly_interactive_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class PlotlyInteractiveReportTests(unittest.TestCase):
    def test_fit_line_returns_expected_summary(self) -> None:
        fit = MODULE.fit_line([0, 1, 2, 3], [1.0, 1.9, 3.1, 4.0])
        self.assertGreater(fit["slope"], 0.9)
        self.assertGreater(fit["r_squared"], 0.98)

    def test_runtime_render_creates_html_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            html_out = Path(tmp_dir) / "report.html"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--html-out",
                    str(html_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertTrue(html_out.exists())
            self.assertEqual(payload["trace_count"], 2)
            self.assertGreater(payload["html_size_bytes"], 1000)


if __name__ == "__main__":
    unittest.main()
