from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "plant-science" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "agriculture-food-and-plant-science"
    / "plantcv-plant-phenotyping-starter"
    / "scripts"
    / "run_plantcv_plant_phenotyping.py"
)


class PlantCvPlantPhenotypingTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_mask_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            image_path = Path(tmp_dir) / "toy.png"
            mask_path = Path(tmp_dir) / "mask.png"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--image-out",
                    str(image_path),
                    "--mask-out",
                    str(mask_path),
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
            self.assertTrue(image_path.exists())
            self.assertTrue(mask_path.exists())
            self.assertEqual(payload["foreground_pixel_count"], 768)
            self.assertEqual(payload["bbox"], [20, 16, 43, 47])
            self.assertEqual(payload["bbox_width"], 24)
            self.assertEqual(payload["bbox_height"], 32)

    def test_rejects_invalid_threshold(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--threshold", "300"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--threshold must be between 0 and 255", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
