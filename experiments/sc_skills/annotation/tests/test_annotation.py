from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"
VALIDATOR = SKILL_DIR / "scripts" / "validate_outputs.py"
TOY_INPUT = SKILL_DIR / "examples" / "toy_input.json"
REFERENCE_EXPORT = SKILL_DIR / "examples" / "reference_metadata_export.tsv"
RUNTIME_CHECK = SKILL_DIR / "scripts" / "check_runtime.py"
ENV_PLANNER = SKILL_DIR / "scripts" / "plan_census_env.py"
REFERENCE_RANKER = SKILL_DIR / "scripts" / "rank_reference_metadata.py"


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_public_workflow_surrogates(self) -> None:
        expectations = json.loads(TOY_INPUT.read_text(encoding="utf-8"))["expected_invariants"]

        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                ["python3", str(SCRIPT), "--outdir", tmpdir],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            validated = subprocess.run(
                ["python3", str(VALIDATOR), "--outdir", tmpdir, "--input", str(TOY_INPUT)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(validated.returncode, 0, validated.stderr)

            outdir = Path(tmpdir)
            annotation = pd.read_csv(outdir / "annotation_table.tsv", sep="\t").set_index("cell_id")
            self.assertEqual(annotation.shape[0], expectations["query_cell_count"])
            self.assertEqual(annotation["predicted_label"].value_counts().to_dict(), expectations["predicted_label_counts"])
            self.assertEqual(annotation["marker_status"].value_counts().to_dict(), expectations["marker_status_counts"])
            self.assertEqual(annotation.loc["cell_001", "predicted_label"], "luminal epithelial")
            self.assertEqual(annotation.loc["cell_003", "predicted_label"], "basal epithelial")
            self.assertEqual(annotation.loc["cell_005", "predicted_label"], "basal epithelial")
            self.assertEqual(annotation.loc["cell_006", "predicted_label"], "fibroblast")
            self.assertEqual(annotation.loc["cell_003", "marker_status"], "borderline")
            self.assertEqual(annotation.loc["cell_005", "marker_status"], "conflict")
            self.assertIn("marker_conflict", annotation.loc["cell_005", "notes"])
            self.assertNotIn("niche", "\t".join(annotation.columns).lower())

            exploration = json.loads((outdir / "qc_dataset_exploration.json").read_text(encoding="utf-8"))
            self.assertEqual(exploration["platform_classification"], "single_cell_spatial")
            self.assertEqual(exploration["selected_batch_column"], expectations["selected_batch_column"])
            self.assertEqual(exploration["query_string"], expectations["query_string"])
            self.assertEqual(exploration["query_variants"][0], expectations["query_string"])
            self.assertEqual(exploration["candidate_batch_columns"][0]["column"], "sample_id")
            self.assertGreater(exploration["candidate_batch_columns"][0]["score"], exploration["candidate_batch_columns"][1]["score"])
            self.assertEqual(exploration["shared_gene_count_after_selection"], expectations["shared_gene_count"])
            self.assertEqual(exploration["reference_only_genes"], expectations["reference_only_genes"])
            self.assertGreater(exploration["explained_variance"][0], 0.6)

            search = pd.read_csv(outdir / "qc_reference_search.tsv", sep="\t")
            self.assertEqual(
                search["dataset_id"].tolist(),
                [
                    "cxg_human_lung_normal_001",
                    "cxg_human_lung_inflamed_002",
                    "cxg_mouse_lung_normal_003",
                ],
            )
            self.assertTrue(bool(search.iloc[0]["selected"]))
            self.assertTrue(bool(search.iloc[0]["accept_for_transfer"]))
            self.assertEqual(search.iloc[0]["rejection_reasons"], "none")
            self.assertIn("disease_not_normal:inflamed", search.iloc[1]["rejection_reasons"])
            self.assertIn("species_mismatch:mouse", search.iloc[2]["rejection_reasons"])
            self.assertGreater(float(search.iloc[0]["score"]), float(search.iloc[1]["score"]))

            selection = json.loads((outdir / "qc_reference_selection.json").read_text(encoding="utf-8"))
            self.assertEqual(selection["selected_dataset_id"], expectations["selected_reference_dataset_id"])
            self.assertEqual(selection["shared_gene_count"], expectations["shared_gene_count"])
            self.assertEqual(selection["reference_only_genes"], expectations["reference_only_genes"])
            self.assertEqual(selection["query_variants"][0], expectations["query_string"])
            self.assertIn("cellxgene_census.download_source_h5ad", selection["download_contract"]["download_call"])
            self.assertEqual(selection["materialized_reference_batches"], {"donor_1": 5, "donor_2": 5})
            self.assertIn("healthy-reference preference", selection["selection_reason"][0])

            label_transfer = pd.read_csv(outdir / "qc_label_transfer.tsv", sep="\t").set_index("cell_id")
            self.assertEqual(label_transfer.loc["cell_001", "broad_label"], "epithelial")
            self.assertEqual(label_transfer.loc["cell_006", "broad_label"], "stromal")
            self.assertTrue(np.all(label_transfer["label_confidence"].to_numpy(dtype=float) > 0.0))
            self.assertAlmostEqual(float(label_transfer.loc["cell_001", "broad_margin"]), 1.0, places=6)
            self.assertGreater(float(label_transfer.loc["cell_006", "broad_margin"]), 0.7)
            self.assertLess(float(label_transfer.loc["cell_005", "subtype_margin"]), 0.55)
            self.assertEqual(label_transfer.loc["cell_003", "marker_status"], "borderline")
            self.assertEqual(label_transfer.loc["cell_005", "marker_status"], "conflict")

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(
                summary["method_steps"],
                [
                    "dataset_exploration_and_platform_detection",
                    "normalize_total_log1p_and_pca_surrogate",
                    "cellxgene_reference_search_and_materialization",
                    "harmony_surrogate_label_transfer",
                    "hierarchical_marker_review",
                ],
            )
            self.assertEqual(summary["selected_reference_dataset_id"], expectations["selected_reference_dataset_id"])
            self.assertEqual(summary["prediction_counts"], expectations["predicted_label_counts"])
            self.assertEqual(summary["marker_status_counts"], expectations["marker_status_counts"])
            self.assertIn("review_flagged_cells", summary)
            self.assertIn("run_summary.json", summary["written_files"])

    def test_spot_based_inputs_are_rejected_and_routed_to_deconvolution(self) -> None:
        payload = json.loads(TOY_INPUT.read_text(encoding="utf-8"))
        payload["query_metadata"]["platform_name"] = "Visium"
        payload["query_metadata"]["observation_unit"] = "spot"
        runner = SKILL_DIR / "scripts" / "run_annotation.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "spot_input.json"
            outdir = Path(tmpdir) / "out"
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            completed = subprocess.run(
                ["python3", str(runner), "--input", str(input_path), "--outdir", str(outdir)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertNotEqual(completed.returncode, 0)
        combined = "\n".join([completed.stdout, completed.stderr])
        self.assertIn("single-cell-resolution spatial data", combined)
        self.assertIn("spatial-deconvolution", combined)

    def test_runtime_check_reports_missing_modules_and_fallbacks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "runtime.json"
            completed = subprocess.run(
                [
                    "python3",
                    str(RUNTIME_CHECK),
                    "--modules",
                    "json",
                    "definitely_missing_annotation_pkg",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(out_path.read_text(encoding="utf-8"))

        checks = {entry["module"]: entry for entry in payload["checked_modules"]}
        self.assertTrue(checks["json"]["available"])
        self.assertFalse(checks["definitely_missing_annotation_pkg"]["available"])
        self.assertIn("pip install definitely_missing_annotation_pkg", checks["definitely_missing_annotation_pkg"]["install_hint"])
        self.assertIn("host_python_support_for_cellxgene_census", payload["fallback_guidance"][0]["missing_if_any"])
        self.assertIn("cellxgene_census_requires_python", payload)
        self.assertIn("host_python_supported_for_cellxgene_census", payload)
        self.assertGreaterEqual(len(payload["bootstrap_commands"]), 2)

    def test_env_planner_routes_python_313_hosts_to_supported_runtime(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "env_plan.json"
            completed = subprocess.run(
                [
                    "python3",
                    str(ENV_PLANNER),
                    "--python-version",
                    "3.13.0",
                    "--conda",
                    "yes",
                    "--python312",
                    "no",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(out_path.read_text(encoding="utf-8"))

        self.assertFalse(payload["host_python_supported_for_cellxgene_census"])
        self.assertEqual(payload["recommended_route"], "conda_py312_census_env")
        self.assertIn("python=3.12", " ".join(payload["routes"]["conda_py312_census_env"]["commands"]))
        self.assertIn("Discover", " ".join(payload["routes"]["discover_export_fallback"]["commands"]))

    def test_reference_ranking_helper_prefers_healthy_matching_annotated_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            ranked_path = Path(tmpdir) / "ranked.tsv"
            summary_path = Path(tmpdir) / "summary.json"
            completed = subprocess.run(
                [
                    "python3",
                    str(REFERENCE_RANKER),
                    "--metadata",
                    str(REFERENCE_EXPORT),
                    "--species",
                    "human",
                    "--tissue",
                    "lung",
                    "--condition",
                    "normal",
                    "--out",
                    str(ranked_path),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            ranked = pd.read_csv(ranked_path, sep="\t")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["selected_dataset_id"], "cxg_human_lung_normal_001")
        self.assertEqual(summary["accepted_candidate_count"], 1)
        self.assertTrue(bool(ranked.iloc[0]["selected"]))
        self.assertTrue(bool(ranked.iloc[0]["accept_for_transfer"]))
        self.assertIn("disease_not_normal:inflamed", ranked.loc[ranked["dataset_id"] == "cxg_human_lung_inflamed_002", "rejection_reasons"].iloc[0])
        self.assertIn("missing_broad_labels", ranked.loc[ranked["dataset_id"] == "cxg_human_lung_sparse_004", "rejection_reasons"].iloc[0])


if __name__ == "__main__":
    unittest.main()
