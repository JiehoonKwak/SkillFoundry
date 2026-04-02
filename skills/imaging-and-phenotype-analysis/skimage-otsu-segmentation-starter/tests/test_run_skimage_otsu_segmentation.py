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
    / "imaging-and-phenotype-analysis"
    / "skimage-otsu-segmentation-starter"
    / "scripts"
    / "run_skimage_otsu_segmentation.py"
)


def load_script_module():
    spec = importlib.util.spec_from_file_location("skimage_otsu_segmentation", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class SkimageOtsuSegmentationTests(unittest.TestCase):
    def test_module_exposes_expected_entry_points(self) -> None:
        self.assertTrue(hasattr(MODULE, "generate_toy_image"))
        self.assertTrue(hasattr(MODULE, "segment_image"))

    def test_runtime_segmentation_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "segmentation.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["image_shape"], [96, 96])
            self.assertEqual(payload["object_count"], 3)
            self.assertGreater(payload["foreground_pixels"], 500)


if __name__ == "__main__":
    unittest.main()
