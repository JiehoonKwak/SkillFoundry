from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CHEMTOOLS_PYTHON = ROOT / "slurm" / "envs" / "chemtools" / "bin" / "python"


class DeepchemFeaturizationTests(unittest.TestCase):
    def test_graph_summary(self) -> None:
        self.assertTrue(CHEMTOOLS_PYTHON.exists(), f"Missing chemtools interpreter: {CHEMTOOLS_PYTHON}")
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "deepchem-summary.json"
            subprocess.run(
                [
                    str(CHEMTOOLS_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/scripts/featurize_molecules.py",
                    "--input",
                    "skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/examples/molecules.tsv",
                    "--out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["molecule_count"], 2)
            self.assertGreater(payload["graphs"][0]["node_count"], 0)
            self.assertGreater(payload["graphs"][0]["node_feature_count"], 0)


if __name__ == "__main__":
    unittest.main()
