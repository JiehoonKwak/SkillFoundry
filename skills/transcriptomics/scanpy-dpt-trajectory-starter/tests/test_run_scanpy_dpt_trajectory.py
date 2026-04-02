from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"
SCRIPT = ROOT / "skills" / "transcriptomics" / "scanpy-dpt-trajectory-starter" / "scripts" / "run_scanpy_dpt_trajectory.py"


class ScanpyDptTrajectoryTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "trajectory.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    "skills/transcriptomics/scanpy-dpt-trajectory-starter/examples/toy_counts.tsv",
                    "--root-cell",
                    "c0",
                    "--expected-order",
                    "skills/transcriptomics/scanpy-dpt-trajectory-starter/examples/expected_order.txt",
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
            self.assertEqual(payload["trajectory_order"], ["c0", "c1", "c2", "c3", "c4", "c5"])
            self.assertTrue(payload["monotonic_in_expected_order"])
            self.assertEqual(payload["pseudotime"][0]["pseudotime"], 0.0)
            self.assertEqual(payload["pseudotime"][-1]["pseudotime"], 1.0)

    def test_rejects_unknown_root_cell(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--counts",
                "skills/transcriptomics/scanpy-dpt-trajectory-starter/examples/toy_counts.tsv",
                "--root-cell",
                "missing-cell",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Root cell", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
