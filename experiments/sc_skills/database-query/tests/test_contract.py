from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_toy_input_stay_aligned(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        expected_deliverables = sorted(item["path"] for item in metadata["deliverables"])
        self.assertEqual(
            expected_deliverables,
            ["query_results.tsv", "resolution_notes.md", "source_provenance.md"],
        )
        self.assertNotIn("deliverables", example)
        self.assertEqual(len(example["queries"]), example["expected_invariants"]["query_count"])
        self.assertEqual(
            sorted(item["query"] for item in example["expected_invariants"]["expected_results"]),
            sorted(example["queries"]),
        )
        self.assertIn("concepts", example["local_lookup_tables"])
        self.assertIn("aliases", example["local_lookup_tables"])
        self.assertIn("mock_api_payloads", example)


if __name__ == "__main__":
    unittest.main()
