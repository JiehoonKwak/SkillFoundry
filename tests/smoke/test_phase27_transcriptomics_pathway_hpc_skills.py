from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCANPY_PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"


class Phase27TranscriptomicsPathwayHpcSmokeTests(unittest.TestCase):
    def test_scanpy_cell_type_annotation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "annotation.json"
            subprocess.run(
                [
                    str(SCANPY_PYTHON),
                    "skills/transcriptomics/scanpy-cell-type-annotation-starter/scripts/run_scanpy_cell_type_annotation.py",
                    "--counts",
                    "skills/transcriptomics/scanpy-cell-type-annotation-starter/examples/toy_counts.tsv",
                    "--markers",
                    "skills/transcriptomics/scanpy-cell-type-annotation-starter/examples/toy_markers.json",
                    "--truth",
                    "skills/transcriptomics/scanpy-cell-type-annotation-starter/examples/toy_truth.tsv",
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
            self.assertEqual(payload["accuracy"], 1.0)
            self.assertEqual(payload["predicted_label_counts"], {"B": 3, "T": 3})

    def test_scanpy_dpt_trajectory_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "trajectory.json"
            subprocess.run(
                [
                    str(SCANPY_PYTHON),
                    "skills/transcriptomics/scanpy-dpt-trajectory-starter/scripts/run_scanpy_dpt_trajectory.py",
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
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["trajectory_order"], ["c0", "c1", "c2", "c3", "c4", "c5"])
            self.assertTrue(payload["monotonic_in_expected_order"])

    def test_reactome_pathway_hierarchy_walk_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "hierarchy.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/reactome-pathway-hierarchy-walk-starter/scripts/run_reactome_hierarchy_walk.py",
                    "--species",
                    "9606",
                    "--stable-id",
                    "R-HSA-141409",
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
            self.assertEqual(payload["top_level_pathway"], "Cell Cycle")
            self.assertEqual(payload["ancestor_count"], 5)

    def test_slurm_job_array_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "array.json"
            subprocess.run(
                [
                    "python3",
                    "skills/hpc/slurm-job-array-starter/scripts/run_slurm_job_array.py",
                    "--partition",
                    "cpu",
                    "--job-name",
                    "phase27-array",
                    "--array-spec",
                    "0-1",
                    "--sleep",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["task_count"], 2)
            self.assertEqual(sorted(task["accounting"]["State"] for task in payload["task_summaries"]), ["COMPLETED", "COMPLETED"])


if __name__ == "__main__":
    unittest.main()
