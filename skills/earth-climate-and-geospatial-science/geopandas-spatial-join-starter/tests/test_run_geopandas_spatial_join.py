from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "geospatial" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "earth-climate-and-geospatial-science"
    / "geopandas-spatial-join-starter"
    / "scripts"
    / "run_geopandas_spatial_join.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("geopandas_spatial_join", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class GeoPandasSpatialJoinTests(unittest.TestCase):
    def test_module_exposes_expected_entry_points(self) -> None:
        self.assertTrue(hasattr(MODULE, "configure_geospatial_env"))
        self.assertTrue(hasattr(MODULE, "summarize_join"))

    def test_runtime_spatial_join_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "spatial-join.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=90,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["assigned_counts"], {"north": 2, "south": 1})
            self.assertEqual(payload["unassigned_point_ids"], ["P4"])
            self.assertEqual(payload["projected_crs"], "EPSG:3857")
            self.assertTrue(str(payload["environment"]["proj_lib"]).endswith("/share/proj"))


if __name__ == "__main__":
    unittest.main()
