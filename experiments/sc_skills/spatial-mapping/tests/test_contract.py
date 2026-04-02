from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_contract_and_toy_input_are_aligned(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        self.assertEqual(
            [item["path"] for item in metadata["deliverables"]],
            ["mapped_labels.tsv", "mapping_scores.tsv", "spatial_mapping_report.md"],
        )
        self.assertEqual(
            metadata["starter_qc_files"],
            [
                "intermediate_qc.json",
                "gene_intersection_qc.tsv",
                "marker_qc.tsv",
                "niche_qc.tsv",
                "run_summary.json",
            ],
        )
        self.assertNotIn("deliverables", example)
        self.assertEqual(len(example["reference"]["profiles"]), 2)
        self.assertEqual(len(example["spatial"]["spots"]), 6)
        self.assertEqual(example["niche"]["n_neighbors"], 3)
        self.assertEqual(example["expected_invariants"]["shared_gene_count"], 8)
        self.assertEqual(
            set(example["expected_invariants"]["primary_labels"]),
            {spot["observation_id"] for spot in example["spatial"]["spots"]},
        )


if __name__ == "__main__":
    unittest.main()
