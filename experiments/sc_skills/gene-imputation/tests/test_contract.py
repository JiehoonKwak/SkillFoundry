from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_toy_input_keeps_raw_inputs_only(self) -> None:
        payload = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", payload)
        self.assertNotIn("mapping_matrix", payload)
        self.assertEqual(len(payload["shared_genes"]), 6)
        self.assertEqual(len(payload["heldout_genes"]), 2)
        self.assertEqual(len(payload["spatial"]["spots"]), 3)
        self.assertEqual(payload["expected_invariants"]["gene_intersection_size"], 6)

    def test_metadata_deliverables_stay_stable(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        observed = [item["path"] for item in metadata["deliverables"]]
        self.assertEqual(
            observed,
            ["imputed_expression.h5ad", "heldout_metrics.tsv", "imputation_report.md"],
        )


if __name__ == "__main__":
    unittest.main()
