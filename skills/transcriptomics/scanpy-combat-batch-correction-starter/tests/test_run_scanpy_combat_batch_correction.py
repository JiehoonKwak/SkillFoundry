from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "scanpy-combat-batch-correction-starter"
    / "scripts"
    / "run_scanpy_combat_batch_correction.py"
)
COUNTS = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "scanpy-combat-batch-correction-starter"
    / "examples"
    / "toy_counts.tsv"
)
METADATA = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "scanpy-combat-batch-correction-starter"
    / "examples"
    / "toy_metadata.tsv"
)


class ScanpyCombatBatchCorrectionTests(unittest.TestCase):
    def test_runtime_summary_shrinks_batch_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "scanpy_combat.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    str(COUNTS),
                    "--metadata",
                    str(METADATA),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["cells"], 8)
            self.assertEqual(payload["genes"], 4)
            self.assertGreater(payload["pre_batch_mean_abs_diff"], payload["post_batch_mean_abs_diff"])
            self.assertLess(payload["batch_gap_ratio"], 0.05)
            self.assertGreater(payload["cell_type_centroid_distance_post"], 1.0)

    def test_rejects_missing_metadata_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_metadata = Path(tmp_dir) / "bad_metadata.tsv"
            bad_metadata.write_text("cell\tbatch\nx\tb1\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    str(COUNTS),
                    "--metadata",
                    str(bad_metadata),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=180,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Missing required columns", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
