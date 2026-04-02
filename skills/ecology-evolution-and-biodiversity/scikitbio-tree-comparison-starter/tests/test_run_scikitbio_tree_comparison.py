from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "ecology" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "ecology-evolution-and-biodiversity"
    / "scikitbio-tree-comparison-starter"
    / "scripts"
    / "run_scikitbio_tree_comparison.py"
)


class ScikitBioTreeComparisonTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_distances(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "tree.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["shared_tip_count"], 4)
            self.assertEqual(payload["left_total_branch_length"], 6.0)
            self.assertEqual(payload["right_total_branch_length"], 6.0)
            self.assertEqual(payload["robinson_foulds_distance"], 4.0)
            self.assertEqual(payload["weighted_robinson_foulds_distance"], 4.0)


if __name__ == "__main__":
    unittest.main()
