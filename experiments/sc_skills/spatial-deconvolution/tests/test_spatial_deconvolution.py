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
    def test_toy_run_computes_deconvolution_metrics_and_qc(self) -> None:
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
            abundance = ad.read_h5ad(outdir / "cell_type_abundance.h5ad")
            selected = np.asarray(abundance.X, dtype=float)
            self.assertEqual(selected.shape, (4, 3))
            self.assertTrue(np.allclose(selected.sum(axis=1), 1.0, atol=1e-6))
            self.assertEqual(abundance.obs_names.tolist(), ["spot_001", "spot_002", "spot_003", "spot_004"])
            self.assertEqual(abundance.var_names.tolist(), ["T_cell", "myeloid", "fibroblast"])
            self.assertEqual(abundance.uns["selected_method"], "cell2location_like")
            self.assertTrue(np.allclose(selected, np.asarray(abundance.layers["cell2location_like"]), atol=1e-6))
            self.assertIn("destvi_like", abundance.layers)
            self.assertIn("uniform_baseline", abundance.layers)
            self.assertIn("spatial", abundance.obsm)

            dominant = abundance.obs["dominant_label"].to_dict()
            self.assertEqual(dominant["spot_001"], "T_cell")
            self.assertEqual(dominant["spot_002"], "myeloid")
            self.assertEqual(dominant["spot_003"], "fibroblast")
            self.assertEqual(dominant["spot_004"], "T_cell")
            self.assertGreater(float(abundance.obs.loc["spot_004", "dominant_score"]), 0.5)

            comparison = pd.read_csv(outdir / "model_comparison.tsv", sep="\t")
            metric_table = (
                comparison.pivot(index="model", columns="metric", values="value")
                .sort_index()
            )
            self.assertLess(
                float(metric_table.loc["cell2location_like", "mean_abundance_rmse"]),
                float(metric_table.loc["destvi_like", "mean_abundance_rmse"]),
            )
            self.assertLess(
                float(metric_table.loc["destvi_like", "mean_abundance_rmse"]),
                float(metric_table.loc["uniform_baseline", "mean_abundance_rmse"]),
            )
            self.assertGreater(
                float(metric_table.loc["cell2location_like", "mean_profile_correlation"]),
                0.99,
            )
            self.assertEqual(float(metric_table.loc["cell2location_like", "selected_output"]), 1.0)

            diagnostics = pd.read_csv(outdir / "qc_spot_diagnostics.tsv", sep="\t")
            self.assertEqual(set(diagnostics["method"]), {"cell2location_like", "destvi_like", "uniform_baseline"})
            self.assertTrue(np.allclose(diagnostics["abundance_sum"], 1.0, atol=1e-6))
            selected_diag = diagnostics.loc[diagnostics["method"] == "cell2location_like"].set_index("spot_id")
            self.assertEqual(selected_diag.loc["spot_004", "dominant_label"], "T_cell")
            self.assertGreater(float(selected_diag.loc["spot_001", "dominant_margin"]), 0.6)
            self.assertGreater(float(selected_diag.loc["spot_002", "dominant_margin"]), 0.4)
            self.assertGreater(float(selected_diag.loc["spot_003", "dominant_margin"]), 0.6)

            neighbors = pd.read_csv(outdir / "qc_neighbor_graph.tsv", sep="\t")
            self.assertEqual(set(neighbors["relation"]), {"self", "neighbor"})
            row_sums = neighbors.groupby("source_spot")["weight"].sum().to_dict()
            self.assertEqual(set(row_sums), {"spot_001", "spot_002", "spot_003", "spot_004"})
            for value in row_sums.values():
                self.assertAlmostEqual(float(value), 1.0, places=6)

            gene_qc = json.loads((outdir / "qc_gene_intersection.json").read_text(encoding="utf-8"))
            self.assertEqual(gene_qc["shared_gene_count"], 5)
            self.assertEqual(gene_qc["cell_types"], ["T_cell", "myeloid", "fibroblast"])

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["selected_method"], "cell2location_like")
            self.assertEqual(
                summary["method_steps"],
                [
                    "gene_intersection",
                    "cell2location_like_nnls",
                    "destvi_like_smoothed_log_nnls",
                    "uniform_baseline_comparison",
                ],
            )
            self.assertEqual(summary["target_spot_dominant_label"], "T_cell")


if __name__ == "__main__":
    unittest.main()
