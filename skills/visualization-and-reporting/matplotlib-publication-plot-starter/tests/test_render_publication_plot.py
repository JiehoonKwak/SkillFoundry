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
    / "matplotlib-publication-plot-starter"
    / "scripts"
    / "render_publication_plot.py"
)
INPUT = (
    ROOT
    / "skills"
    / "visualization-and-reporting"
    / "matplotlib-publication-plot-starter"
    / "examples"
    / "toy_measurements.tsv"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("render_publication_plot", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class MatplotlibPublicationPlotTests(unittest.TestCase):
    def test_fit_line_returns_expected_summary(self) -> None:
        fit = MODULE.fit_line([0, 1, 2, 3], [1.0, 2.0, 2.8, 4.1])
        self.assertGreater(fit["slope"], 0.9)
        self.assertLess(fit["slope"], 1.1)
        self.assertGreater(fit["r_squared"], 0.95)

    def test_format_x_tick_preserves_non_integer_labels(self) -> None:
        self.assertEqual(MODULE.format_x_tick(1.5), "1.5")
        self.assertEqual(MODULE.format_x_tick(2.0), "2")

    def test_runtime_render_creates_png_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            png_out = Path(tmp_dir) / "figure.png"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--png-out",
                    str(png_out),
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
            self.assertTrue(png_out.exists())
            self.assertEqual(payload["point_count"], 6)
            self.assertGreater(payload["r_squared"], 0.98)


if __name__ == "__main__":
    unittest.main()
