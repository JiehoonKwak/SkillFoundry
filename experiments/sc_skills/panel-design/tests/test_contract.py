from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_toy_input_contains_raw_inputs_not_deliverables(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        toy_input = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        deliverable_paths = sorted(item["path"] for item in metadata["deliverables"])
        self.assertEqual(
            deliverable_paths,
            ["panel_candidates.tsv", "panel_rationale.md", "platform_notes.md"],
        )
        self.assertNotIn("deliverables", toy_input)
        self.assertEqual(len(toy_input["candidate_markers"]), toy_input["expected_invariants"]["aggregated_candidate_count"])
        self.assertEqual(len(toy_input["cell_types"]), 5)
        self.assertEqual(
            sum(int(item["quota"]) for item in toy_input["cell_types"]),
            5,
        )
        self.assertGreaterEqual(toy_input["panel_size"], 6)


if __name__ == "__main__":
    unittest.main()
