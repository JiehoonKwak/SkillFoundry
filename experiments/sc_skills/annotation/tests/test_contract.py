from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_toy_input_stays_raw_and_aligned_with_contract(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        example = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        expected = example["expected_invariants"]

        self.assertEqual(
            sorted(item["path"] for item in metadata["deliverables"]),
            ["annotation_table.tsv", "marker_evidence.md"],
        )
        self.assertEqual(
            metadata["starter_qc_files"],
            [
                "qc_dataset_exploration.json",
                "qc_reference_search.tsv",
                "qc_reference_selection.json",
                "qc_label_transfer.tsv",
                "run_summary.json",
            ],
        )
        self.assertEqual(
            [item["path"] for item in metadata["helper_scripts"]],
            [
                "scripts/check_runtime.py",
                "scripts/plan_census_env.py",
                "scripts/rank_reference_metadata.py",
                "scripts/run_annotation.py",
            ],
        )
        self.assertIn("00b_census_env_plan.json", metadata["real_run_artifacts"])
        self.assertIn("04_reference_ranked.tsv", metadata["real_run_artifacts"])
        self.assertIn("05_reference_selection.json", metadata["real_run_artifacts"])

        for forbidden_key in ["annotation_table.tsv", "marker_evidence.md", "deliverables"]:
            self.assertNotIn(forbidden_key, example)

        self.assertEqual(len(example["query"]["cells"]), expected["query_cell_count"])
        self.assertEqual(len(example["toy_czi_catalog"]), expected["catalog_candidate_count"])
        self.assertEqual(len(example["query"]["genes"]), expected["shared_gene_count"])
        self.assertEqual(len(example["marker_panels"]["subtype"]), expected["subtype_count"])
        self.assertEqual(len({entry["label"] for entry in example["marker_panels"]["broad"]}), expected["broad_class_count"])
        self.assertEqual(example["reference_selection"]["min_cells"], 1000)

        catalog_ids = {item["dataset_id"] for item in example["toy_czi_catalog"]}
        self.assertIn(expected["selected_reference_dataset_id"], catalog_ids)
        selected_row = next(item for item in example["toy_czi_catalog"] if item["dataset_id"] == expected["selected_reference_dataset_id"])
        self.assertIn(selected_row["reference_key"], example["reference_datasets"])

        selected_reference = example["reference_datasets"][selected_row["reference_key"]]
        self.assertEqual(selected_reference["dataset_id"], expected["selected_reference_dataset_id"])
        self.assertEqual(
            sorted(set(selected_reference["genes"]) - set(example["query"]["genes"])),
            expected["reference_only_genes"],
        )
        self.assertIn(expected["selected_batch_column"], example["query"]["cells"][0])

        for cell in example["query"]["cells"]:
            self.assertEqual(len(cell["counts"]), len(example["query"]["genes"]))
        for cell in selected_reference["cells"]:
            self.assertEqual(len(cell["counts"]), len(selected_reference["genes"]))

    def test_runbook_and_reference_export_schema_cover_real_run_path(self) -> None:
        skill_text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        assets_text = (SKILL_DIR / "assets" / "README.md").read_text(encoding="utf-8")

        for required_phrase in [
            "Step 0. Check the runtime before touching the data",
            "Step 0b. Choose the Census/Discover route explicitly",
            "Step 1. Inspect the query dataset and classify the platform",
            "scripts/check_runtime.py",
            "scripts/plan_census_env.py",
            "scripts/rank_reference_metadata.py",
            "cellxgene_census.download_source_h5ad",
            "scanpy.external.pp.harmony_integrate",
            "scanpy.tl.ingest",
            "spatial-deconvolution",
            "tissue-niche-annotation",
            "Python 3.13",
            "Deterministic toy validation",
        ]:
            self.assertIn(required_phrase, skill_text)

        self.assertIn("Real-run path", assets_text)
        self.assertIn("Deterministic surrogate path", assets_text)
        self.assertIn("conda", assets_text)
        self.assertIn("Discover browser", assets_text)

        reference_export = pd.read_csv(SKILL_DIR / "examples" / "reference_metadata_export.tsv", sep="\t")
        self.assertEqual(
            reference_export.columns.tolist(),
            [
                "dataset_id",
                "dataset_title",
                "collection_name",
                "species",
                "tissue",
                "disease",
                "assay",
                "suspension_type",
                "dataset_total_cell_count",
                "has_broad_labels",
                "has_subtype_labels",
                "cell_type_obs_key",
                "subtype_obs_key",
                "is_primary_data",
                "dataset_h5ad_uri",
            ],
        )
        self.assertGreaterEqual(reference_export.shape[0], 4)
        self.assertIn("cxg_human_lung_normal_001", reference_export["dataset_id"].tolist())


if __name__ == "__main__":
    unittest.main()
