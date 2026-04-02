from __future__ import annotations

import json
import unittest
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_raw_toy_input_stay_aligned(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        self.assertEqual(
            sorted(item["path"] for item in metadata["deliverables"]),
            ["neighbor_enrichment.tsv", "spatial_stats.tsv", "squidpy_report.md"],
        )
        self.assertEqual(
            metadata["starter_qc_files"],
            [
                "qc_spatial_graph.tsv",
                "qc_label_scores.tsv",
                "qc_pairwise_stats.tsv",
                "run_summary.json",
            ],
        )
        self.assertNotIn("deliverables", example)
        self.assertEqual(example["expected_invariants"]["spot_count"], len(example["spots"]))
        self.assertEqual(example["expected_invariants"]["gene_count"], len(example["genes"]))
        self.assertEqual(example["expected_invariants"]["image_feature_count"], len(example["image_features"]))
        self.assertEqual(example["expected_invariants"]["domain_count"], len(example["candidate_domains"]))

        for domain in example["candidate_domains"]:
            self.assertEqual(sorted(example["domain_signatures"][domain].keys()), sorted(example["genes"]))
            self.assertEqual(sorted(example["domain_image_profiles"][domain].keys()), sorted(example["image_features"]))

        for spot in example["spots"]:
            self.assertEqual(sorted(spot["counts"].keys()), sorted(example["genes"]))
            self.assertEqual(sorted(spot["image_features"].keys()), sorted(example["image_features"]))
            self.assertNotIn("label", spot)


if __name__ == "__main__":
    unittest.main()
