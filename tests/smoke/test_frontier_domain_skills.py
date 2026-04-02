from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
MATERIALS_PYTHON = ROOT / "slurm" / "envs" / "materials" / "bin" / "python"
GEOSPATIAL_PYTHON = ROOT / "slurm" / "envs" / "geospatial" / "bin" / "python"
STATISTICS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class FrontierDomainSkillSmokeTests(unittest.TestCase):
    def test_skimage_otsu_segmentation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "skimage.json"
            subprocess.run(
                [
                    str(SCIENTIFIC_PYTHON),
                    "skills/imaging-and-phenotype-analysis/skimage-otsu-segmentation-starter/scripts/run_skimage_otsu_segmentation.py",
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

    def test_xarray_climate_cube_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "xarray.json"
            subprocess.run(
                [
                    str(SCIENTIFIC_PYTHON),
                    "skills/earth-climate-and-geospatial-science/xarray-climate-cube-starter/scripts/run_xarray_climate_cube.py",
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
            self.assertEqual(payload["dims"], {"lat": 3, "lon": 4, "time": 12})
            self.assertEqual(payload["hottest_month"], "2025-04")

    def test_matminer_composition_featurization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "matminer.json"
            subprocess.run(
                [
                    str(MATERIALS_PYTHON),
                    "skills/materials-science-and-engineering/matminer-composition-featurization/scripts/run_matminer_composition_features.py",
                    "--formula",
                    "Fe2O3",
                    "--formula",
                    "LiFePO4",
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
            self.assertEqual(payload["formula_count"], 2)
            self.assertEqual(payload["summaries"][0]["reduced_formula"], "Fe2O3")

    def test_geopandas_spatial_join_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "geopandas.json"
            subprocess.run(
                [
                    str(GEOSPATIAL_PYTHON),
                    "skills/earth-climate-and-geospatial-science/geopandas-spatial-join-starter/scripts/run_geopandas_spatial_join.py",
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
            self.assertEqual(payload["assigned_counts"], {"north": 2, "south": 1})
            self.assertEqual(payload["unassigned_point_ids"], ["P4"])

    def test_pymc_bayesian_linear_regression_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "pymc.json"
            subprocess.run(
                [
                    str(STATISTICS_PYTHON),
                    "skills/statistical-and-machine-learning-foundations-for-science/pymc-bayesian-linear-regression-starter/scripts/run_pymc_linear_regression.py",
                    "--input",
                    "skills/statistical-and-machine-learning-foundations-for-science/pymc-bayesian-linear-regression-starter/examples/toy_observations.tsv",
                    "--draws",
                    "120",
                    "--tune",
                    "120",
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
            self.assertGreater(payload["posterior_means"]["slope"], 0.7)
            self.assertLess(payload["posterior_means"]["slope"], 1.0)


if __name__ == "__main__":
    unittest.main()
