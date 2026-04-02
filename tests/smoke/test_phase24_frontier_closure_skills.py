from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
METAGENOMICS_PYTHON = ROOT / "slurm" / "envs" / "metagenomics" / "bin" / "python"
GEOSPATIAL_PYTHON = ROOT / "slurm" / "envs" / "geospatial" / "bin" / "python"
CAUSAL_PYTHON = ROOT / "slurm" / "envs" / "causal" / "bin" / "python"
TRANSCRIPTOMICS_PYTHON = ROOT / "slurm" / "envs" / "transcriptomics" / "bin" / "python"


class Phase24FrontierClosureSmokeTests(unittest.TestCase):
    def test_sourmash_signature_compare_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "sourmash.json"
            subprocess.run(
                [
                    str(METAGENOMICS_PYTHON),
                    "skills/genomics/sourmash-signature-compare-starter/scripts/run_sourmash_signature_compare.py",
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
            self.assertEqual(payload["shared_hash_count"], 12)
            self.assertAlmostEqual(payload["jaccard_similarity"], 0.631579, places=6)

    def test_rasterio_windowed_raster_preprocessing_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "rasterio.json"
            subprocess.run(
                [
                    str(GEOSPATIAL_PYTHON),
                    "skills/earth-climate-and-geospatial-science/rasterio-windowed-raster-preprocessing-starter/scripts/run_rasterio_windowed_preprocessing.py",
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
            self.assertEqual(payload["window_values"], [[6.0, 7.0], [10.0, 11.0]])

    def test_dowhy_average_treatment_effect_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "dowhy.json"
            subprocess.run(
                [
                    str(CAUSAL_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/dowhy-average-treatment-effect-starter/scripts/run_dowhy_average_treatment_effect.py",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertAlmostEqual(payload["average_treatment_effect"], 1.855648, places=6)
            self.assertLess(abs(payload["placebo_effect"]), 0.1)

    def test_pydeseq2_differential_expression_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "pydeseq2.json"
            subprocess.run(
                [
                    str(TRANSCRIPTOMICS_PYTHON),
                    "skills/transcriptomics/pydeseq2-differential-expression-starter/scripts/run_pydeseq2_differential_expression.py",
                    "--counts",
                    "skills/transcriptomics/pydeseq2-differential-expression-starter/examples/toy_counts.tsv",
                    "--metadata",
                    "skills/transcriptomics/pydeseq2-differential-expression-starter/examples/toy_metadata.tsv",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["top_gene"], "GeneD")
            self.assertEqual(payload["significant_gene_count"], 4)

    def test_quarto_notebook_report_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            html_out = Path(tmp_dir) / "report.html"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/visualization-and-reporting/quarto-notebook-report-starter/scripts/render_quarto_notebook_report.py",
                    "--input",
                    "skills/visualization-and-reporting/quarto-notebook-report-starter/examples/toy_report.ipynb",
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
            self.assertTrue(payload["title_present"])
            self.assertTrue(payload["executed_output_present"])
            self.assertTrue(html_out.exists())


if __name__ == "__main__":
    unittest.main()
