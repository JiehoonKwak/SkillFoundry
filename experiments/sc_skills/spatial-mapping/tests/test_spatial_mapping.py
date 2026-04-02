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


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_mapping_qc_and_review_invariants(self) -> None:
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

            expectations = json.loads(TOY_INPUT.read_text(encoding="utf-8"))["expected_invariants"]
            outdir = Path(tmpdir)

            mapped = pd.read_csv(outdir / "mapped_labels.tsv", sep="\t").set_index("observation_id")
            self.assertEqual(mapped["primary_label"].to_dict(), expectations["primary_labels"])
            self.assertEqual(mapped["marker_review_status"].to_dict(), expectations["review_statuses"])
            self.assertTrue(np.allclose(mapped["gene_overlap_count"], expectations["shared_gene_count"]))
            self.assertTrue(np.allclose(mapped["gene_overlap_fraction"], 0.8))
            self.assertGreater(mapped.loc["spot_001", "score_margin"], 0.08)
            self.assertLess(
                mapped.loc["spot_003", "score_margin"],
                expectations["ambiguous_observations"]["spot_003"]["max_score_margin"],
            )
            self.assertLess(
                mapped.loc["spot_006", "score_margin"],
                expectations["ambiguous_observations"]["spot_006"]["max_score_margin"],
            )

            scores = pd.read_csv(outdir / "mapping_scores.tsv", sep="\t")
            global_scores = scores.loc[scores["scope"] == "global"].set_index("metric_name")
            self.assertGreaterEqual(
                float(global_scores.loc["shared_gene_fraction", "value"]),
                expectations["metrics"]["shared_gene_fraction_min"],
            )
            self.assertGreaterEqual(
                float(global_scores.loc["mean_primary_score", "value"]),
                expectations["metrics"]["mean_primary_score_min"],
            )
            self.assertGreaterEqual(
                float(global_scores.loc["mean_score_margin", "value"]),
                expectations["metrics"]["mean_score_margin_min"],
            )
            self.assertGreaterEqual(
                float(global_scores.loc["niche_agreement_rate", "value"]),
                expectations["metrics"]["niche_agreement_rate_min"],
            )
            self.assertGreaterEqual(
                float(global_scores.loc["supported_fraction", "value"]),
                expectations["metrics"]["supported_fraction_min"],
            )

            observation_scores = scores.loc[scores["scope"] == "observation"]
            margin_status = observation_scores.loc[observation_scores["metric_name"] == "score_margin"].set_index("target_id")
            self.assertEqual(margin_status.loc["spot_003", "qc_status"], "review")
            self.assertEqual(margin_status.loc["spot_006", "qc_status"], "review")

            gene_qc = pd.read_csv(outdir / "gene_intersection_qc.tsv", sep="\t")
            self.assertEqual(gene_qc.loc[gene_qc["role"] == "shared"].shape[0], expectations["shared_gene_count"])
            self.assertEqual(
                sorted(gene_qc.loc[gene_qc["role"] == "reference_only", "gene"].tolist()),
                expectations["reference_only_genes"],
            )
            self.assertEqual(
                sorted(gene_qc.loc[gene_qc["role"] == "spatial_only", "gene"].tolist()),
                expectations["spatial_only_genes"],
            )

            marker_qc = pd.read_csv(outdir / "marker_qc.tsv", sep="\t").set_index("observation_id")
            self.assertGreater(marker_qc.loc["spot_001", "marker_gap"], 1.5)
            self.assertGreater(
                marker_qc.loc["spot_003", "marker_gap"],
                expectations["ambiguous_observations"]["spot_003"]["min_marker_gap"],
            )
            self.assertLess(
                marker_qc.loc["spot_006", "marker_gap"],
                expectations["ambiguous_observations"]["spot_006"]["max_marker_gap"],
            )

            niche_qc = pd.read_csv(outdir / "niche_qc.tsv", sep="\t").set_index("observation_id")
            self.assertTrue((niche_qc["niche_agreement"] == 1.0).all())
            self.assertEqual(niche_qc.loc["spot_003", "niche_majority_label"], "epithelial")
            self.assertEqual(niche_qc.loc["spot_006", "niche_majority_label"], "immune")

            intermediate = json.loads((outdir / "intermediate_qc.json").read_text(encoding="utf-8"))
            similarity = np.asarray(intermediate["cosine_similarity"]["matrix"], dtype=float)
            self.assertEqual(similarity.shape, (6, 2))
            self.assertGreater(float(similarity[0].max()), 0.99)
            self.assertLess(float(abs(similarity[2, 0] - similarity[2, 1])), 0.03)

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(
                summary["method_steps"],
                [
                    "gene_intersection",
                    "normalize_total_log1p",
                    "reference_centroid_cosine_similarity",
                    "top_label_margin",
                    "3nn_niche_majority_vote",
                    "marker_consistency_review",
                ],
            )
            self.assertEqual(summary["prediction_counts"], {"epithelial": 3, "immune": 3})
            self.assertEqual(summary["review_counts"], {"mixed": 1, "needs_review": 1, "supported": 4})


if __name__ == "__main__":
    unittest.main()
