from __future__ import annotations

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
    / "rasterio-windowed-raster-preprocessing-starter"
    / "scripts"
    / "run_rasterio_windowed_preprocessing.py"
)


class RasterioWindowedPreprocessingTests(unittest.TestCase):
    def test_runtime_summary_contains_expected_window(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["input_shape"], [4, 4])
            self.assertEqual(payload["window_values"], [[6.0, 7.0], [10.0, 11.0]])
            self.assertEqual(payload["resampled_values"], [[3.5, 5.5], [11.5, 13.5]])


if __name__ == "__main__":
    unittest.main()
