from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_pseudotime_velocity_and_fate_invariants(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                ["python3", str(SCRIPT), "--outdir", tmpdir],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            outdir = Path(tmpdir)
            example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
            trajectory = pd.read_csv(outdir / "trajectory_coordinates.tsv", sep="\t").set_index("cell_id")

            self.assertEqual(len(trajectory), example["expected_invariants"]["cell_count"])
            self.assertEqual(trajectory.loc["cell_00", "lineage_label"], "root")
            self.assertEqual(
                trajectory["lineage_label"].value_counts().to_dict(),
                {"root": 1, "trunk": 3, "terminal_alpha": 3, "terminal_beta": 3},
            )
            self.assertAlmostEqual(float(trajectory.loc["cell_00", "pseudotime"]), 0.0, places=6)
            self.assertAlmostEqual(float(trajectory.loc["cell_06", "pseudotime"]), 1.0, places=6)
            self.assertAlmostEqual(float(trajectory.loc["cell_09", "pseudotime"]), 1.0, places=6)

            for path_key in ("alpha_path", "beta_path"):
                path = example["expected_invariants"][path_key]
                path_times = [float(trajectory.loc[cell_id, "pseudotime"]) for cell_id in path]
                self.assertTrue(all(left < right for left, right in zip(path_times, path_times[1:])))

            fate = pd.read_csv(outdir / "fate_probabilities.tsv", sep="\t").pivot(
                index="cell_id",
                columns="terminal_state",
                values="probability",
            )
            self.assertTrue(np.allclose(fate.sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6))
            self.assertGreater(float(fate.loc["cell_05", "terminal_alpha"]), 0.99)
            self.assertGreater(float(fate.loc["cell_08", "terminal_beta"]), 0.99)
            self.assertEqual(float(fate.loc["cell_06", "terminal_alpha"]), 1.0)
            self.assertEqual(float(fate.loc["cell_09", "terminal_beta"]), 1.0)
            self.assertAlmostEqual(
                float(fate.loc["cell_03", "terminal_alpha"]),
                float(fate.loc["cell_03", "terminal_beta"]),
                delta=0.05,
            )

            velocity = pd.read_csv(outdir / "qc_velocity_surrogate.tsv", sep="\t").set_index("cell_id")
            self.assertGreater(float(velocity.loc["cell_00", "velocity_dx"]), 0.0)
            self.assertGreater(float(velocity.loc["cell_04", "velocity_dy"]), 0.0)
            self.assertLess(float(velocity.loc["cell_07", "velocity_dy"]), 0.0)
            self.assertEqual(int(velocity.loc["cell_06", "forward_neighbor_count"]), 0)
            self.assertEqual(int(velocity.loc["cell_09", "forward_neighbor_count"]), 0)

            transition = pd.read_csv(outdir / "qc_transition_matrix.tsv", sep="\t")
            transition_sums = transition.groupby("source_cell")["transition_probability"].sum().to_numpy(dtype=float)
            self.assertTrue(np.allclose(transition_sums, 1.0, atol=1e-6))

            preprocessing = json.loads((outdir / "qc_preprocessing.json").read_text(encoding="utf-8"))
            self.assertTrue(preprocessing["connected"])
            self.assertEqual(preprocessing["component_count"], 1)
            self.assertEqual(preprocessing["root_cell_id"], "cell_00")

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["branchpoint_cell_id"], "cell_03")
            self.assertEqual(
                summary["method_steps"],
                [
                    "knn_graph_construction",
                    "shortest_path_pseudotime",
                    "neighbor_delta_velocity_surrogate",
                    "absorbing_markov_fate_probabilities",
                ],
            )


if __name__ == "__main__":
    unittest.main()
