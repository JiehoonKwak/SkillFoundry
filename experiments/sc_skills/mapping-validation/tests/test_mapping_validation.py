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
TOY_INPUT = SKILL_DIR / "examples" / "toy_input.json"


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_expected_invariants(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                ["python3", str(SCRIPT), "--outdir", tmpdir],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            expectations = json.loads(TOY_INPUT.read_text(encoding="utf-8"))["expected_invariants"]
            metrics = pd.read_csv(Path(tmpdir) / "validation_metrics.tsv", sep="\t")
            metric_map = dict(zip(metrics["metric"], metrics["value"]))

            self.assertGreaterEqual(metric_map["heldout_pearson_r"], expectations["metrics"]["heldout_pearson_r_min"])
            self.assertLessEqual(metric_map["heldout_rmse"], expectations["metrics"]["heldout_rmse_max"])
            self.assertGreaterEqual(metric_map["cv_mean_pearson_r"], expectations["metrics"]["cv_mean_pearson_r_min"])
            self.assertLessEqual(metric_map["cv_mean_rmse"], expectations["metrics"]["cv_mean_rmse_max"])
            self.assertGreaterEqual(metric_map["marker_recovery_rate"], expectations["metrics"]["marker_recovery_rate_min"])
            self.assertLessEqual(
                metric_map["mean_assignment_entropy"],
                expectations["metrics"]["mean_assignment_entropy_max"],
            )
            self.assertGreaterEqual(
                metric_map["neighbor_label_agreement_rate"],
                expectations["metrics"]["neighbor_label_agreement_rate_min"],
            )
            self.assertGreaterEqual(
                metric_map["mean_top_cosine_similarity"],
                expectations["metrics"]["mean_top_cosine_similarity_min"],
            )

            qc = json.loads((Path(tmpdir) / "intermediate_qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["gene_intersection"]["count"], expectations["gene_intersection_size"])
            self.assertEqual(
                {item["spot_id"]: item["top_cell_type"] for item in qc["spot_assignments"]},
                expectations["top_assignments"],
            )

            weight_rows = np.asarray(qc["mapping_weights"], dtype=float)
            self.assertEqual(weight_rows.shape, (3, 2))
            self.assertTrue(np.allclose(weight_rows.sum(axis=1), 1.0, atol=1e-6))
            self.assertEqual(len(qc["cross_validation"]), 6)
            self.assertTrue(all(item["marker_recovered"] for item in qc["marker_recovery"]))

            gene_cv = pd.read_csv(Path(tmpdir) / "gene_cv_qc.tsv", sep="\t")
            self.assertEqual(gene_cv.shape[0], 6)
            self.assertTrue((gene_cv["pearson_r"] >= 0.99).all())

            assignments = pd.read_csv(Path(tmpdir) / "spot_assignment_qc.tsv", sep="\t")
            self.assertEqual(list(assignments["top_cell_type"]), ["epithelial", "epithelial", "immune"])
            self.assertTrue((assignments["dominant_margin"] > 0).all())

            coherence = pd.read_csv(Path(tmpdir) / "spatial_coherence_qc.tsv", sep="\t")
            self.assertEqual(coherence.shape[0], 3)
            self.assertTrue(coherence["agreement_with_neighbors"].between(0.0, 1.0).all())

            holdout_text = (Path(tmpdir) / "holdout_summary.md").read_text(encoding="utf-8")
            self.assertIn("## Held-out genes", holdout_text)
            report_text = (Path(tmpdir) / "mapping_validation_report.md").read_text(encoding="utf-8")
            self.assertIn("## Validation metrics", report_text)


if __name__ == "__main__":
    unittest.main()
