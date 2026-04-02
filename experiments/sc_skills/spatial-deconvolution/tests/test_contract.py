from __future__ import annotations

import json
import unittest
from pathlib import Path


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_toy_input_define_raw_starter_inputs(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))

        self.assertEqual(
            [item["path"] for item in metadata["deliverables"]],
            ["cell_type_abundance.h5ad", "model_comparison.tsv", "deconvolution_summary.md"],
        )
        self.assertNotIn("deliverables", example)

        reference_genes = example["reference"]["genes"]
        spatial_genes = example["spatial"]["genes"]
        shared_genes = [gene for gene in spatial_genes if gene in set(reference_genes)]
        self.assertEqual(len(shared_genes), example["expected_invariants"]["shared_gene_count"])

        cell_types = [profile["label"] for profile in example["reference"]["profiles"]]
        self.assertEqual(len(cell_types), example["expected_invariants"]["cell_type_count"])
        self.assertEqual(example["synthetic_truth"]["cell_types"], cell_types)

        for profile in example["reference"]["profiles"]:
            self.assertEqual(len(profile["centroid"]), len(reference_genes))

        expected_spot_count = example["expected_invariants"]["spot_count"]
        self.assertEqual(len(example["spatial"]["spots"]), expected_spot_count)
        self.assertEqual(len(example["synthetic_truth"]["abundances"]), expected_spot_count)

        truth_by_spot = {row["spot_id"]: row["weights"] for row in example["synthetic_truth"]["abundances"]}
        self.assertEqual(set(truth_by_spot), set(example["expected_invariants"]["dominant_labels"]))
        for spot in example["spatial"]["spots"]:
            self.assertEqual(len(spot["counts"]), len(spatial_genes))
            self.assertGreater(sum(spot["counts"]), 0)
            weights = truth_by_spot[spot["spot_id"]]
            self.assertAlmostEqual(sum(weights), 1.0, places=6)
            dominant_index = max(range(len(weights)), key=weights.__getitem__)
            self.assertEqual(
                cell_types[dominant_index],
                example["expected_invariants"]["dominant_labels"][spot["spot_id"]],
            )


if __name__ == "__main__":
    unittest.main()
