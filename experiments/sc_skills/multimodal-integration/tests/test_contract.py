from __future__ import annotations

import json
import unittest
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_keeps_declared_deliverables(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        expected = [item["path"] for item in metadata["deliverables"]]
        self.assertEqual(
            expected,
            ["integrated_latent.h5mu", "modality_qc.tsv", "integration_report.md"],
        )

    def test_toy_input_contains_raw_modalities_not_precomputed_deliverables(self) -> None:
        payload = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", payload)
        self.assertEqual(len(payload["cells"]), 8)
        for cell in payload["cells"]:
            self.assertEqual(len(cell["rna"]), len(payload["rna_features"]))
            self.assertEqual(len(cell["protein"]), len(payload["protein_features"]))
            self.assertEqual(len(cell["chromatin"]), len(payload["chromatin_features"]))
        self.assertEqual(set(payload["expected_invariants"]["query_labels"]), {"cell_005", "cell_006", "cell_007", "cell_008"})


if __name__ == "__main__":
    unittest.main()
