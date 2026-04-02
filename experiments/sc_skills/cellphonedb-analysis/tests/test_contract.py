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
                "cellphonedb_significant_means.tsv",
                "interaction_ranked.tsv",
                "cellphonedb_report.md",
            ],
        )
        observed_qc = sorted(item["path"] for item in metadata["qc_artifacts"])
        self.assertEqual(
            observed_qc,
            sorted(
                [
                    "qc_group_expression.tsv",
                    "qc_statistical_null.tsv",
                    "qc_deg_support.tsv",
                    "qc_adjacency_support.tsv",
                    "run_summary.json",
                ]
            ),
        )

    def test_toy_input_contains_raw_inputs_not_precomputed_outputs(self) -> None:
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", example)
        self.assertNotIn("cellphonedb_report.md", example)
        self.assertEqual(len(example["cells"]), example["expected_invariants"]["cell_count"])
        self.assertEqual(len(example["ligand_receptor_pairs"]), example["expected_invariants"]["pair_count"])
        self.assertEqual(len({cell["group"] for cell in example["cells"]}), example["expected_invariants"]["group_count"])
        self.assertGreater(len(example["degs"]), 0)
        self.assertGreater(len(example["adjacency"]), 0)


if __name__ == "__main__":
    unittest.main()
