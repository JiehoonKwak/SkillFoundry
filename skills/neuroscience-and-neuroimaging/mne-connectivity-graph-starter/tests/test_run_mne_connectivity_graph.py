from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "neuro" / "bin" / "python"
SCRIPT = ROOT / "skills" / "neuroscience-and-neuroimaging" / "mne-connectivity-graph-starter" / "scripts" / "run_mne_connectivity_graph.py"


class MneConnectivityGraphTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_edge_structure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "connectivity_graph.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["channels"], ["Fz", "Cz", "Pz"])
            self.assertEqual(payload["epochs"], 8)
            self.assertEqual(payload["max_edge"], {"source": "Fz", "target": "Cz", "weight": 0.926824})
            self.assertEqual(payload["strong_edges"], [{"source": "Fz", "target": "Cz", "weight": 0.926824}])
            self.assertEqual(payload["connectivity_matrix"][2][0], 0.119457)


if __name__ == "__main__":
    unittest.main()
