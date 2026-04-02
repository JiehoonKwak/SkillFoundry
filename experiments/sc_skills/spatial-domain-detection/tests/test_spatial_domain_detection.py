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


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_graph_derived_domains_and_markers(self) -> None:
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
            labels = pd.read_csv(outdir / "domain_labels.tsv", sep="\t").set_index("cell_id")
            self.assertEqual(labels.shape[0], 8)
            self.assertEqual(labels["selected_method"].nunique(), 1)
            self.assertEqual(labels["selected_method"].iat[0], "spagcn_like")
            self.assertEqual(labels["domain_label"].value_counts().to_dict(), {"epithelial": 3, "immune": 3, "stromal": 2})
            self.assertEqual(labels.loc["spot_003", "base_best_domain"], "stromal")
            self.assertEqual(labels.loc["spot_003", "domain_label"], "epithelial")
            self.assertEqual(labels.loc["spot_006", "base_best_domain"], "stromal")
            self.assertEqual(labels.loc["spot_006", "domain_label"], "immune")
            self.assertGreater(float(labels.loc["spot_003", "domain_score"]), 0.6)
            self.assertGreater(float(labels.loc["spot_006", "confidence_margin"]), 0.08)

            markers = pd.read_csv(outdir / "domain_markers.tsv", sep="\t")
            marker_lookup = {
                (row.domain_label, row.marker_gene): float(row.effect_size)
                for row in markers.itertuples()
            }
            self.assertGreater(marker_lookup[("epithelial", "EPCAM")], 1.3)
            self.assertGreater(marker_lookup[("stromal", "COL1A1")], 1.1)
            self.assertGreater(marker_lookup[("immune", "CXCL13")], 1.4)

            domain_scores = pd.read_csv(outdir / "qc_domain_scores.tsv", sep="\t")
            self.assertEqual(set(domain_scores["method"]), {"spagcn_like", "graphst_like"})
            self.assertEqual(domain_scores.groupby("method")["cell_id"].nunique().to_dict(), {"graphst_like": 8, "spagcn_like": 8})
            spagcn_scores = domain_scores.loc[domain_scores["method"] == "spagcn_like"].set_index("cell_id")
            graphst_scores = domain_scores.loc[domain_scores["method"] == "graphst_like"].set_index("cell_id")
            self.assertTrue(bool(spagcn_scores.loc["spot_003", "graph_shifted"]))
            self.assertTrue(bool(spagcn_scores.loc["spot_006", "graph_shifted"]))
            self.assertEqual(graphst_scores.loc["spot_006", "assigned_domain"], "stromal")
            self.assertGreater(float(graphst_scores.loc["spot_006", "domain_score"]), 0.95)

            graph = pd.read_csv(outdir / "qc_spatial_graph.tsv", sep="\t")
            self.assertEqual(set(graph["relation"]), {"neighbor", "self"})
            row_sums = graph.groupby(["method", "source_cell"])["weight"].sum().to_numpy(dtype=float)
            self.assertTrue(np.allclose(row_sums, 1.0, atol=1e-6))

            comparison = pd.read_csv(outdir / "qc_method_comparison.tsv", sep="\t")
            metric_table = comparison.pivot(index="method", columns="metric", values="value")
            self.assertGreater(
                float(metric_table.loc["spagcn_like", "composite_score"]),
                float(metric_table.loc["graphst_like", "composite_score"]),
            )
            self.assertGreater(
                float(metric_table.loc["spagcn_like", "weighted_edge_agreement"]),
                float(metric_table.loc["graphst_like", "weighted_edge_agreement"]),
            )
            self.assertEqual(float(metric_table.loc["spagcn_like", "selected_output"]), 1.0)
            self.assertEqual(float(metric_table.loc["graphst_like", "selected_output"]), 0.0)

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["selected_method"], "spagcn_like")
            self.assertEqual(summary["boundary_spots"], ["spot_003", "spot_006"])
            self.assertEqual(
                summary["method_steps"],
                [
                    "spatial_graph_construction",
                    "signature_scoring",
                    "spagcn_like_refinement",
                    "graphst_like_partitioning",
                ],
            )


if __name__ == "__main__":
    unittest.main()
