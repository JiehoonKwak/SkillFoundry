from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "synthetic-toy-dataset-generator-starter"
SCRIPT = SKILL / "scripts" / "generate_synthetic_toy_dataset.py"


class SyntheticToyDatasetTests(unittest.TestCase):
    def run_skill(self, *extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *extra],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_generates_deterministic_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle_dir = Path(tmpdir) / "bundle"
            summary_path = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--sample-count", "6",
                "--feature-count", "4",
                "--seed", "17",
                "--out-dir", str(bundle_dir),
                "--summary-out", str(summary_path),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["matrix_shape"], [4, 6])
            self.assertEqual(payload["matrix_checksum"], 1306)
            self.assertTrue((bundle_dir / "samples.tsv").exists())
            self.assertTrue((bundle_dir / "features.tsv").exists())
            self.assertTrue((bundle_dir / "matrix.tsv").exists())
            self.assertTrue((bundle_dir / "manifest.json").exists())

    def test_same_seed_repeats_same_checksum(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle_a = Path(tmpdir) / "bundle_a"
            bundle_b = Path(tmpdir) / "bundle_b"
            summary_a = Path(tmpdir) / "summary_a.json"
            summary_b = Path(tmpdir) / "summary_b.json"
            args = [
                "--sample-count", "5",
                "--feature-count", "3",
                "--seed", "42",
            ]
            result_a = self.run_skill(*args, "--out-dir", str(bundle_a), "--summary-out", str(summary_a))
            result_b = self.run_skill(*args, "--out-dir", str(bundle_b), "--summary-out", str(summary_b))
            self.assertEqual(result_a.returncode, 0, result_a.stderr)
            self.assertEqual(result_b.returncode, 0, result_b.stderr)
            payload_a = json.loads(summary_a.read_text(encoding="utf-8"))
            payload_b = json.loads(summary_b.read_text(encoding="utf-8"))
            self.assertEqual(payload_a["matrix_checksum"], payload_b["matrix_checksum"])
            self.assertEqual(payload_a["samples_preview"], payload_b["samples_preview"])

    def test_rejects_invalid_sample_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle_dir = Path(tmpdir) / "bundle"
            summary_path = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--sample-count", "1",
                "--feature-count", "4",
                "--seed", "17",
                "--out-dir", str(bundle_dir),
                "--summary-out", str(summary_path),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("sample_count", result.stderr)


if __name__ == "__main__":
    unittest.main()
