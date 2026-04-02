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
            ["domain_detection_report.md", "domain_labels.tsv", "domain_markers.tsv"],
        )
        self.assertEqual(
            metadata["starter_qc_files"],
            [
                "qc_spatial_graph.tsv",
                "qc_domain_scores.tsv",
                "qc_method_comparison.tsv",
                "run_summary.json",
            ],
        )
        self.assertNotIn("deliverables", example)
        self.assertEqual(example["expected_invariants"]["spot_count"], len(example["spots"]))
        self.assertEqual(example["expected_invariants"]["gene_count"], len(example["genes"]))
        self.assertEqual(example["expected_invariants"]["domain_count"], len(example["candidate_domains"]))

        for domain in example["candidate_domains"]:
            self.assertEqual(sorted(example["domain_signatures"][domain].keys()), sorted(example["genes"]))
        for spot in example["spots"]:
            self.assertEqual(sorted(spot["counts"].keys()), sorted(example["genes"]))
            self.assertEqual(len(spot["histology_rgb"]), 3)


if __name__ == "__main__":
    unittest.main()
