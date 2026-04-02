from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHEM_TOOLS_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
SCIENTIFIC_PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
STATISTICS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
MATERIALS_PYTHON = ROOT / "slurm" / "envs" / "materials" / "bin" / "python"


class Phase23FrontierSkillSmokeTests(unittest.TestCase):
    def test_skimage_regionprops_feature_extraction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "regionprops.json"
            subprocess.run(
                [
                    str(SCIENTIFIC_PYTHON),
                    "skills/imaging-and-phenotype-analysis/skimage-regionprops-feature-extraction/scripts/run_skimage_regionprops_features.py",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["object_count"], 3)

    def test_rdkit_molecule_standardization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "rdkit-standardized.json"
            subprocess.run(
                [
                    str(CHEM_TOOLS_PYTHON),
                    "skills/drug-discovery-and-cheminformatics/rdkit-molecule-standardization/scripts/standardize_rdkit_molecule.py",
                    "--smiles",
                    "CC(=O)[O-].[Na+]",
                    "--name",
                    "sodium-acetate",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=90,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["canonical_tautomer_smiles"], "CC(=O)O")

    def test_arviz_posterior_diagnostics_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "arviz.json"
            subprocess.run(
                [
                    str(STATISTICS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/arviz-posterior-diagnostics-starter/scripts/run_arviz_posterior_diagnostics.py",
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
            self.assertLess(payload["max_rhat"], 1.02)

    def test_plotly_interactive_report_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            html_out = Path(tmp_dir) / "report.html"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(STATISTICS_PYTHON),
                    "skills/visualization-and-reporting/plotly-interactive-report-starter/scripts/render_plotly_interactive_report.py",
                    "--input",
                    "skills/visualization-and-reporting/plotly-interactive-report-starter/examples/toy_measurements.tsv",
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
            self.assertEqual(payload["trace_count"], 2)
            self.assertTrue(html_out.exists())

    def test_precommit_regression_testing_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workspace = Path(tmp_dir) / "workspace"
            out_path = Path(tmp_dir) / "precommit.json"
            subprocess.run(
                [
                    "python3",
                    "skills/meta-maintenance/precommit-regression-testing-starter/scripts/run_precommit_regression.py",
                    "--workspace",
                    str(workspace),
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
            self.assertEqual(payload["returncode"], 0)

    def test_matminer_property_regression_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "matminer-regression.json"
            subprocess.run(
                [
                    str(MATERIALS_PYTHON),
                    "skills/materials-science-and-engineering/matminer-property-regression-starter/scripts/run_matminer_property_regression.py",
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
            self.assertEqual(payload["test_count"], 2)
            self.assertGreater(payload["feature_count"], 100)


if __name__ == "__main__":
    unittest.main()
