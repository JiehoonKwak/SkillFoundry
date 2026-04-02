from __future__ import annotations

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
    / "skimage-regionprops-feature-extraction"
    / "scripts"
    / "run_skimage_regionprops_features.py"
)


class SkimageRegionpropsFeatureExtractionTests(unittest.TestCase):
    def test_runtime_summary_contains_three_objects(self) -> None:
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
            self.assertEqual(payload["object_count"], 3)
            self.assertEqual(payload["threshold"], 0.35)
            self.assertGreater(payload["objects"][0]["area"], payload["objects"][-1]["area"])
            self.assertGreater(payload["mean_intensity_range"][1], payload["mean_intensity_range"][0])

    def test_rejects_non_positive_size(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--size", "0"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("size must be a positive integer", completed.stderr + completed.stdout)

    def test_rejects_invalid_threshold(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--threshold", "1.2"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("threshold must be between 0 and 1", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
