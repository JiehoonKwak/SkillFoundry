from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "systems-biology" / "networkx-network-propagation-starter" / "scripts" / "run_networkx_network_propagation.py"
INPUT = ROOT / "skills" / "systems-biology" / "networkx-network-propagation-starter" / "examples" / "toy_network.tsv"
SEEDS = ROOT / "skills" / "systems-biology" / "networkx-network-propagation-starter" / "examples" / "toy_seeds.txt"


class NetworkXNetworkPropagationTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_ranking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "network_propagation.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--seeds",
                    str(SEEDS),
                    "--top-k",
                    "5",
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
            self.assertEqual(payload["node_count"], 14)
            self.assertEqual(payload["edge_count"], 13)
            self.assertEqual(payload["seed_nodes"], ["EGFR", "ERBB2"])
            self.assertEqual(payload["top_non_seed_nodes"][0]["node"], "GRB2")
            self.assertAlmostEqual(payload["score_sum"], 1.0, places=5)

    def test_rejects_missing_seed_file(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                str(INPUT),
                "--seeds",
                "missing.txt",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Seed file not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
