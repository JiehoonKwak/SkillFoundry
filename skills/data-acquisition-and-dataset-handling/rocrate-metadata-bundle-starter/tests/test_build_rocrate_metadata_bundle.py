from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "rocrate-metadata-bundle-starter"
    / "scripts"
    / "build_rocrate_metadata_bundle.py"
)
INPUT = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "rocrate-metadata-bundle-starter"
    / "examples"
    / "toy_measurements.csv"
)


class ROCrateMetadataBundleTests(unittest.TestCase):
    def test_builds_crate_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            crate_dir = tmp_root / "toy_bundle"
            out_path = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--crate-dir",
                    str(crate_dir),
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
            self.assertTrue(payload["has_metadata_file"])
            self.assertEqual(payload["root_dataset_name"], "Toy Measurement Bundle")
            self.assertEqual(payload["data_entity_count"], 1)
            self.assertEqual(payload["bundled_files"], ["toy_measurements.csv"])

    def test_missing_input_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(Path(tmp_dir) / "missing.csv"),
                    "--crate-dir",
                    str(Path(tmp_dir) / "out"),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input file not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
