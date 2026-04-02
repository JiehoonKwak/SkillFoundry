from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "systems-biology"
    / "networkx-graph-construction-starter"
    / "scripts"
    / "run_networkx_graph_construction.py"
)
INPUT = (
    ROOT
    / "skills"
    / "systems-biology"
    / "networkx-graph-construction-starter"
    / "examples"
    / "toy_pathway_edges.tsv"
)


class NetworkXGraphConstructionTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_topology(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "networkx_graph.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--source-node",
                    "EGFR",
                    "--target-node",
                    "STAT3",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["node_count"], 12)
            self.assertEqual(payload["edge_count"], 11)
            self.assertEqual(payload["connected_component_count"], 1)
            self.assertEqual(payload["top_degree_centrality"][0]["node"], "MAPK1")
            self.assertEqual(
                payload["shortest_path_nodes"],
                ["EGFR", "GRB2", "SOS1", "KRAS", "BRAF", "MAP2K1", "MAPK1", "STAT3"],
            )

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                "missing.tsv",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input edge list not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
