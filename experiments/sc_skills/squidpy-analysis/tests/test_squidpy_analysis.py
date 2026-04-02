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
    def test_toy_run_computes_graph_derived_spatial_statistics(self) -> None:
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
            labels = pd.read_csv(outdir / "qc_label_scores.tsv", sep="\t").set_index("cell_id")
            self.assertEqual(labels.shape[0], 8)
            self.assertEqual(labels["assigned_label"].value_counts().to_dict(), {"epithelial": 3, "immune": 3, "stromal": 2})
            self.assertEqual(sorted(labels.index[labels["graph_shifted"]].tolist()), ["spot_003", "spot_006"])
            self.assertEqual(labels.loc["spot_003", "base_best_domain"], "stromal")
            self.assertEqual(labels.loc["spot_003", "assigned_label"], "epithelial")
            self.assertEqual(labels.loc["spot_006", "base_best_domain"], "stromal")
            self.assertEqual(labels.loc["spot_006", "assigned_label"], "immune")
            self.assertGreater(float(labels.loc["spot_003", "domain_score"]), 0.45)
            self.assertGreater(float(labels.loc["spot_006", "domain_score"]), 0.43)
            self.assertGreater(float(labels.loc["spot_003", "confidence_margin"]), 0.08)
            self.assertGreater(float(labels.loc["spot_006", "confidence_margin"]), 0.05)

            pairwise = pd.read_csv(outdir / "qc_pairwise_stats.tsv", sep="\t").set_index(["label_a", "label_b"])
            self.assertEqual(pairwise.shape[0], 9)
            self.assertGreater(float(pairwise.loc[("stromal", "stromal"), "enrichment_score"]), 1.3)
            self.assertGreater(float(pairwise.loc[("epithelial", "epithelial"), "enrichment_score"]), 1.2)
            self.assertLess(float(pairwise.loc[("epithelial", "immune"), "enrichment_score"]), -4.0)
            self.assertAlmostEqual(float(pairwise.loc[("epithelial", "epithelial"), "co_occurrence_score"]), 2.153846, places=5)
            self.assertAlmostEqual(float(pairwise.loc[("immune", "immune"), "co_occurrence_score"]), 2.153846, places=5)
            self.assertEqual(float(pairwise.loc[("epithelial", "immune"), "co_occurrence_score"]), 0.0)

            stats = pd.read_csv(outdir / "spatial_stats.tsv", sep="\t")
            co_occurrence = stats.loc[stats["statistic"] == "co_occurrence"].set_index("scope")
            moran = stats.loc[stats["statistic"] == "moran_i"].set_index("scope")
            self.assertGreater(float(co_occurrence.loc["epithelial|epithelial", "value"]), float(co_occurrence.loc["epithelial|stromal", "value"]))
            self.assertEqual(float(co_occurrence.loc["epithelial|immune", "value"]), 0.0)
            self.assertGreater(float(moran.loc["gene:EPCAM", "value"]), 0.77)
            self.assertGreater(float(moran.loc["image:epithelial_texture", "value"]), 0.84)
            self.assertGreater(
                float(moran.loc["image:epithelial_texture", "value"]),
                float(moran.loc["image:fibrousness", "value"]),
            )

            graph = pd.read_csv(outdir / "qc_spatial_graph.tsv", sep="\t")
            self.assertEqual(set(graph["relation"]), {"neighbor", "self"})
            row_sums = graph.groupby("source_cell")["weight"].sum().to_numpy(dtype=float)
            self.assertTrue(np.allclose(row_sums, 1.0, atol=1e-6))
            self.assertTrue((graph.groupby("source_cell").size() == 3).all())

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(
                summary["method_steps"],
                ["spatial_neighbors", "nhood_enrichment", "co_occurrence", "morans_i"],
            )
            self.assertEqual(summary["boundary_spots"], ["spot_003", "spot_006"])
            self.assertEqual(summary["derived_domain_counts"], {"epithelial": 3, "immune": 3, "stromal": 2})
            self.assertEqual(summary["top_neighbor_enrichment"], {"label_a": "stromal", "label_b": "stromal", "enrichment_score": 1.34354})
            self.assertEqual(summary["top_moran"], {"scope": "image:epithelial_texture", "value": 0.846876})

            report = (outdir / "squidpy_report.md").read_text(encoding="utf-8")
            self.assertIn("## Run context", report)
            self.assertIn("## Spatial graph", report)
            self.assertIn("## Statistics", report)
            self.assertIn("## Caveats", report)
            self.assertIn("permutation z-score", report)


if __name__ == "__main__":
    unittest.main()
