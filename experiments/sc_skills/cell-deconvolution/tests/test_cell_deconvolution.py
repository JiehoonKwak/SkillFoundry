from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_mapping_projection_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                ["python3", str(SCRIPT), "--outdir", tmpdir],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            outdir = Path(tmpdir)
            mapping = ad.read_h5ad(outdir / "mapping_matrix.h5ad")
            weights = np.asarray(mapping.X, dtype=float)
            self.assertEqual(weights.shape, (3, 2))
            self.assertTrue(np.allclose(weights.sum(axis=1), 1.0, atol=1e-6))
            self.assertEqual(mapping.obs_names.tolist(), ["spot_001", "spot_002", "spot_003"])
            self.assertEqual(mapping.var_names.tolist(), ["T_cell", "myeloid"])

            dominant = mapping.obs["dominant_label"].to_dict()
            self.assertEqual(dominant["spot_001"], "T_cell")
            self.assertEqual(dominant["spot_002"], "myeloid")
            self.assertEqual(dominant["spot_003"], "T_cell")

            entropy = mapping.obs["assignment_entropy"].to_dict()
            self.assertLess(entropy["spot_001"], 0.9)
            self.assertLess(entropy["spot_002"], 0.9)
            self.assertGreater(entropy["spot_003"], 0.99)

            margin = mapping.obs["dominant_margin"].to_dict()
            self.assertGreater(margin["spot_001"], 0.45)
            self.assertGreater(margin["spot_002"], 0.40)
            self.assertLess(margin["spot_003"], 0.1)

            assignments = pd.read_csv(outdir / "cell_assignment.tsv", sep="\t").set_index("spot_id")
            self.assertEqual(assignments.loc["spot_001", "dominant_label"], "T_cell")
            self.assertEqual(assignments.loc["spot_002", "dominant_label"], "myeloid")
            self.assertGreater(assignments.loc["spot_001", "dominant_score"], 0.7)
            self.assertGreater(assignments.loc["spot_002", "dominant_score"], 0.7)
            self.assertLess(assignments.loc["spot_003", "dominant_score"], 0.55)

            gene_qc = json.loads((outdir / "qc_gene_intersection.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_qc["shared_gene_count"], 8)
            self.assertEqual(gene_qc["fit_gene_count"], 6)
            self.assertEqual(gene_qc["holdout_gene_count"], 2)
            self.assertEqual(gene_qc["holdout_genes"], ["LTB", "CTSB"])

            segment_projection = pd.read_csv(outdir / "qc_segment_projection.tsv", sep="\t")
            assigned_counts = segment_projection.groupby("spot_id")["assigned_segment_count"].sum().to_dict()
            self.assertEqual(assigned_counts, {"spot_001": 3, "spot_002": 3, "spot_003": 2})
            segment_by_label = (
                segment_projection.pivot(index="spot_id", columns="label", values="assigned_segment_count")
                .fillna(0)
                .astype(int)
            )
            self.assertEqual(segment_by_label.loc["spot_001", "T_cell"], 2)
            self.assertEqual(segment_by_label.loc["spot_002", "myeloid"], 2)
            self.assertEqual(segment_by_label.loc["spot_003", "T_cell"], 1)
            self.assertEqual(segment_by_label.loc["spot_003", "myeloid"], 1)

            holdout = pd.read_csv(outdir / "qc_holdout_predictions.tsv", sep="\t")
            holdout_metrics = {}
            for gene, subset in holdout.groupby("gene"):
                observed = subset["observed_value"].to_numpy(dtype=float)
                predicted = subset["predicted_value"].to_numpy(dtype=float)
                holdout_metrics[gene] = float(np.corrcoef(observed, predicted)[0, 1])
                self.assertLess(float(subset["absolute_error"].mean()), 0.25)
            self.assertGreater(holdout_metrics["LTB"], 0.99)
            self.assertGreater(holdout_metrics["CTSB"], 0.99)

            spot_metrics = pd.read_csv(outdir / "qc_spot_metrics.tsv", sep="\t").set_index("spot_id")
            self.assertGreater(spot_metrics.loc["spot_001", "dominant_marker_recovery"], 0.7)
            self.assertGreater(spot_metrics.loc["spot_002", "dominant_marker_recovery"], 0.65)
            self.assertLess(spot_metrics.loc["spot_003", "dominant_marker_recovery"], 0.6)

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(
                summary["method_steps"],
                [
                    "gene_intersection",
                    "nnls_mapping",
                    "segmentation_projection",
                    "holdout_gene_validation",
                ],
            )
            self.assertEqual(summary["prediction_counts"], {"T_cell": 2, "myeloid": 1})


if __name__ == "__main__":
    unittest.main()
