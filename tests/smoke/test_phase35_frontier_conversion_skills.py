from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
REPORTING_PYTHON = ROOT / "slurm" / "envs" / "reporting" / "bin" / "python"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class Phase35FrontierConversionSmokeTests(unittest.TestCase):
    def test_zarr_chunked_array_store_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            summary_out = tmp_root / "summary.json"
            store_out = tmp_root / "toy_matrix.zarr"
            subprocess.run(
                [
                    "python3",
                    "skills/data-acquisition-and-dataset-handling/zarr-chunked-array-store-starter/scripts/run_zarr_chunked_array_store.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/zarr-chunked-array-store-starter/examples/toy_matrix.tsv",
                    "--store-out",
                    str(store_out),
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
            self.assertEqual(payload["chunk_shape"], [2, 2])
            self.assertEqual(payload["matrix_sum"], 300.0)

    def test_openmm_langevin_dynamics_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/openmm-langevin-dynamics-starter/scripts/run_openmm_langevin_dynamics.py",
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
            self.assertEqual(payload["platform"], "Reference")
            self.assertEqual(payload["frame_count"], 6)

    def test_optuna_bayesian_optimization_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(STATS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/optuna-bayesian-optimization-starter/scripts/run_optuna_bayesian_optimization.py",
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
            self.assertEqual(payload["trial_count"], 32)
            self.assertLess(payload["best_value"], 0.05)

    def test_papermill_parameterized_notebook_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            summary_out = tmp_root / "summary.json"
            output_notebook = tmp_root / "executed.ipynb"
            subprocess.run(
                [
                    str(REPORTING_PYTHON),
                    "skills/reproducible-workflows/papermill-parameterized-notebook-starter/scripts/run_papermill_parameterized_notebook.py",
                    "--input",
                    "skills/reproducible-workflows/papermill-parameterized-notebook-starter/examples/toy_parameters.ipynb",
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
            self.assertTrue(payload["injected_parameters_present"])
            self.assertEqual(payload["result"]["sum"], 12)

    def test_github_actions_scientific_ci_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            summary_out = tmp_root / "summary.json"
            workflow_out = tmp_root / "sciskill_ci.yml"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/github-actions-scientific-ci-starter/scripts/render_github_actions_scientific_ci.py",
                    "--workflow-out",
                    str(workflow_out),
                    "--summary-out",
                    str(summary_out),
                    "--smoke-target",
                    "smoke-zarr",
                    "--smoke-target",
                    "smoke-openmm-md",
                    "--smoke-target",
                    "smoke-optuna",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertTrue(payload["contains_validate"])
            self.assertEqual(payload["smoke_target_count"], 3)
            self.assertIn("validate-test", payload["job_names"])


if __name__ == "__main__":
    unittest.main()
