from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_deliverables_and_qc_are_declared(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        observed_deliverables = [item["path"] for item in metadata["deliverables"]]
        self.assertEqual(
            observed_deliverables,
            [
                "liana_rankings.tsv",
                "resource_overlap.tsv",
                "liana_report.md",
            ],
        )
        observed_qc = sorted(item["path"] for item in metadata["qc_artifacts"])
        self.assertEqual(
            observed_qc,
            sorted(
                [
                    "qc_metadata_harmonization.tsv",
                    "qc_group_expression.tsv",
                    "qc_method_scores.tsv",
                    "qc_condition_summary.tsv",
                    "qc_adjacency_support.tsv",
                    "run_summary.json",
                ]
            ),
        )

    def test_toy_input_contains_raw_inputs_not_precomputed_outputs(self) -> None:
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", example)
        self.assertNotIn("liana_report.md", example)
        self.assertEqual(len(example["cells"]), example["expected_invariants"]["cell_count"])
        self.assertEqual(len(example["ligand_receptor_catalog"]), example["expected_invariants"]["pair_count"])
        self.assertEqual(len(example["selected_methods"]), example["expected_invariants"]["selected_method_count"])
        self.assertGreater(len(example["adjacency"]), 0)


if __name__ == "__main__":
    unittest.main()
