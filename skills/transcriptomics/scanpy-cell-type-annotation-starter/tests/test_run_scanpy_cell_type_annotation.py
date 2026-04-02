from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"
SCRIPT = ROOT / "skills" / "transcriptomics" / "scanpy-cell-type-annotation-starter" / "scripts" / "run_scanpy_cell_type_annotation.py"


class ScanpyCellTypeAnnotationTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "annotation.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
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
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["cells"], 6)
            self.assertEqual(payload["genes"], 6)
            self.assertEqual(payload["accuracy"], 1.0)
            self.assertEqual(payload["predicted_label_counts"], {"B": 3, "T": 3})
            self.assertEqual(payload["predictions"][0]["predicted_label"], "T")

    def test_rejects_marker_sets_without_matching_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            markers_path = Path(tmp_dir) / "markers.json"
            markers_path.write_text('{"Unknown": ["NOT_A_GENE"]}\n', encoding="utf-8")
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    "skills/transcriptomics/scanpy-cell-type-annotation-starter/examples/toy_counts.tsv",
                    "--markers",
                    str(markers_path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("No marker genes", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
