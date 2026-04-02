from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_deliverables_match_expected_contract(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        expected = [
            "holdout_summary.md",
            "mapping_validation_report.md",
            "validation_metrics.tsv",
        ]
        observed = sorted(item["path"] for item in metadata["deliverables"])
        self.assertEqual(expected, observed)

    def test_toy_input_stays_raw_and_not_precomputed_outputs(self) -> None:
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", example)
        self.assertIn("reference", example)
        self.assertIn("spatial", example)
        self.assertIn("shared_genes", example)
        self.assertIn("heldout_genes", example)
        self.assertTrue(set(example["shared_genes"]).isdisjoint(example["heldout_genes"]))


if __name__ == "__main__":
    unittest.main()
