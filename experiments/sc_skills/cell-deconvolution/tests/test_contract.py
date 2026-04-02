from __future__ import annotations

import json
import unittest
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_toy_input_define_a_runnable_raw_input_contract(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        self.assertEqual(
            [item["path"] for item in metadata["deliverables"]],
            ["mapping_matrix.h5ad", "cell_assignment.tsv", "deconvolution_report.md"],
        )
        self.assertNotIn("deliverables", example)

        shared_genes = [
            gene for gene in example["spatial"]["genes"] if gene in set(example["reference"]["genes"])
        ]
        self.assertEqual(len(shared_genes), 8)
        self.assertEqual(len(example["fit_genes"]), 6)
        self.assertEqual(len(example["holdout_genes"]), 2)
        self.assertFalse(set(example["fit_genes"]) & set(example["holdout_genes"]))
        self.assertEqual(set(shared_genes), set(example["fit_genes"]) | set(example["holdout_genes"]))

        for profile in example["reference"]["profiles"]:
            self.assertEqual(len(profile["centroid"]), len(example["reference"]["genes"]))
        for spot in example["spatial"]["spots"]:
            self.assertEqual(len(spot["counts"]), len(example["spatial"]["genes"]))
            self.assertGreater(spot["segmentation_cell_count"], 0)


if __name__ == "__main__":
    unittest.main()
