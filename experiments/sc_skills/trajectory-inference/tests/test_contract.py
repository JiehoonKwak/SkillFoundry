from __future__ import annotations

import json
import unittest
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_preserves_declared_deliverables(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        self.assertEqual(
            [item["path"] for item in metadata["deliverables"]],
            [
                "trajectory_coordinates.tsv",
                "fate_probabilities.tsv",
                "trajectory_report.md",
            ],
        )
        self.assertEqual(
            metadata["qc_artifacts"],
            [
                "qc_preprocessing.json",
                "qc_knn_graph.tsv",
                "qc_velocity_surrogate.tsv",
                "qc_transition_matrix.tsv",
                "run_summary.json",
            ],
        )

    def test_toy_input_is_raw_only_and_matches_expected_shape(self) -> None:
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", example)
        self.assertIn("cells", example)
        self.assertIn("terminal_states", example)
        self.assertIn("expected_invariants", example)
        self.assertEqual(len(example["cells"]), example["expected_invariants"]["cell_count"])
        self.assertEqual(example["root_cell_id"], example["expected_invariants"]["root_cell_id"])
        self.assertEqual(example["program_names"], ["alpha_program", "beta_program"])
        observed_cells = {cell["cell_id"] for cell in example["cells"]}
        terminal_cells = {state["cell_id"] for state in example["terminal_states"]}
        self.assertTrue(terminal_cells.issubset(observed_cells))


if __name__ == "__main__":
    unittest.main()
