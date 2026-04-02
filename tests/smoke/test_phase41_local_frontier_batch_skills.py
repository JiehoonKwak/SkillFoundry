from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA_PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
MAINT_PYTHON = ROOT / "slurm" / "envs" / "maintenance" / "bin" / "python"
REPORT_PYTHON = ROOT / "slurm" / "envs" / "reporting" / "bin" / "python"


class Phase41LocalFrontierBatchSkillTests(unittest.TestCase):
    def test_numcodecs_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "numcodecs.json"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    "skills/data-acquisition-and-dataset-handling/numcodecs-compression-decompression-starter/scripts/run_numcodecs_compression_decompression.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/numcodecs-compression-decompression-starter/examples/toy_matrix.tsv",
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
            self.assertTrue(payload["decoded_equal"])

    def test_pyarrow_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            parquet_out = Path(tmp_dir) / "toy.parquet"
            out_path = Path(tmp_dir) / "pyarrow.json"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    "skills/data-acquisition-and-dataset-handling/pyarrow-format-conversion-starter/scripts/run_pyarrow_format_conversion.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/pyarrow-format-conversion-starter/examples/toy_matrix.tsv",
                    "--parquet-out",
                    str(parquet_out),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertTrue(payload["roundtrip_rows_equal"])

    def test_statistics_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            scipy_out = Path(tmp_dir) / "scipy.json"
            umap_out = Path(tmp_dir) / "umap.json"
            pydoe_out = Path(tmp_dir) / "pydoe.json"
            subprocess.run(
                [
                    str(STATS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/scipy-statistical-testing-starter/scripts/run_scipy_statistical_testing.py",
                    "--input",
                    "skills/statistical-and-machine-learning-foundations-for-science/scipy-statistical-testing-starter/examples/toy_groups.tsv",
                    "--out",
                    str(scipy_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    str(STATS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/umap-dimensionality-reduction-starter/scripts/run_umap_dimensionality_reduction.py",
                    "--input",
                    "skills/statistical-and-machine-learning-foundations-for-science/umap-dimensionality-reduction-starter/examples/toy_embedding_input.tsv",
                    "--out",
                    str(umap_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    str(STATS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/pydoe3-experimental-design-starter/scripts/run_pydoe3_experimental_design.py",
                    "--input",
                    "skills/statistical-and-machine-learning-foundations-for-science/pydoe3-experimental-design-starter/examples/toy_factors.json",
                    "--out",
                    str(pydoe_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertEqual(json.loads(scipy_out.read_text(encoding="utf-8"))["higher_mean_group"], "treated")
            self.assertEqual(json.loads(umap_out.read_text(encoding="utf-8"))["sample_count"], 6)
            self.assertEqual(json.loads(pydoe_out.read_text(encoding="utf-8"))["row_count"], 4)

    def test_ase_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "ase.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/ase-geometry-optimization-starter/scripts/run_ase_geometry_optimization.py",
                    "--input",
                    "skills/computational-chemistry-and-molecular-simulation/ase-geometry-optimization-starter/examples/toy_argon_dimer.json",
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
            self.assertLess(payload["final_distance"], payload["initial_distance"])

    def test_maintenance_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            rapidfuzz_out = Path(tmp_dir) / "rapidfuzz.json"
            datasketch_out = Path(tmp_dir) / "datasketch.json"
            subprocess.run(
                [
                    str(MAINT_PYTHON),
                    "skills/meta-maintenance/rapidfuzz-skill-deduplication-starter/scripts/run_rapidfuzz_skill_deduplication.py",
                    "--input",
                    "skills/meta-maintenance/rapidfuzz-skill-deduplication-starter/examples/toy_skills.tsv",
                    "--threshold",
                    "85",
                    "--out",
                    str(rapidfuzz_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    str(MAINT_PYTHON),
                    "skills/meta-maintenance/datasketch-resource-deduplication-starter/scripts/run_datasketch_resource_deduplication.py",
                    "--input",
                    "skills/meta-maintenance/datasketch-resource-deduplication-starter/examples/toy_resources.jsonl",
                    "--threshold",
                    "0.5",
                    "--out",
                    str(datasketch_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertGreaterEqual(len(json.loads(rapidfuzz_out.read_text(encoding="utf-8"))["candidate_pairs"]), 2)
            self.assertEqual(len(json.loads(datasketch_out.read_text(encoding="utf-8"))["candidate_pairs"]), 1)

    def test_mkdocs_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workspace = Path(tmp_dir) / "workspace"
            out_path = Path(tmp_dir) / "mkdocs.json"
            subprocess.run(
                [
                    str(REPORT_PYTHON),
                    "skills/visualization-and-reporting/mkdocs-summary-catalog-starter/scripts/build_mkdocs_summary_catalog.py",
                    "--input",
                    "skills/visualization-and-reporting/mkdocs-summary-catalog-starter/examples/toy_catalog.json",
                    "--workspace",
                    str(workspace),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertIn("index.html", payload["html_files"])


if __name__ == "__main__":
    unittest.main()
