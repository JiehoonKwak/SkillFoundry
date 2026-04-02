from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "earth-climate-and-geospatial-science"
    / "xarray-climate-cube-starter"
    / "scripts"
    / "run_xarray_climate_cube.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("xarray_climate_cube", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class XarrayClimateCubeTests(unittest.TestCase):
    def test_module_exposes_expected_entry_points(self) -> None:
        self.assertTrue(hasattr(MODULE, "build_toy_cube"))
        self.assertTrue(hasattr(MODULE, "summarize_cube"))

    def test_runtime_climate_cube_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "climate-cube.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["dims"], {"lat": 3, "lon": 4, "time": 12})
            self.assertEqual(payload["hottest_month"], "2025-04")
            self.assertGreater(payload["seasonal_range_c"], 10.0)


if __name__ == "__main__":
    unittest.main()
