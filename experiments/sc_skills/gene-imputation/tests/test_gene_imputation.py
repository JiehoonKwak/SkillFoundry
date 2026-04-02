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
            metrics = pd.read_csv(Path(tmpdir) / "heldout_metrics.tsv", sep="\t")
            metric_map = dict(zip(metrics["metric"], metrics["value"]))

            self.assertGreaterEqual(metric_map["heldout_pearson_r"], expectations["metrics"]["heldout_pearson_r_min"])
            self.assertLessEqual(metric_map["heldout_rmse"], expectations["metrics"]["heldout_rmse_max"])
            self.assertGreaterEqual(metric_map["marker_recovery_rate"], expectations["metrics"]["marker_recovery_rate_min"])
            self.assertLessEqual(metric_map["mean_assignment_entropy"], expectations["metrics"]["mean_assignment_entropy_max"])

            qc = json.loads((Path(tmpdir) / "intermediate_qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["gene_intersection"]["count"], expectations["gene_intersection_size"])
            self.assertEqual(
                {item["spot_id"]: item["top_cell_type"] for item in qc["spot_assignments"]},
                expectations["top_assignments"],
            )
            weight_rows = np.asarray(qc["mapping_weights"], dtype=float)
            self.assertEqual(weight_rows.shape, (3, 2))
            self.assertTrue(np.allclose(weight_rows.sum(axis=1), 1.0, atol=1e-6))
            self.assertTrue(all(item["marker_recovered"] for item in qc["marker_recovery"]))

            adata = ad.read_h5ad(Path(tmpdir) / "imputed_expression.h5ad")
            self.assertEqual(adata.n_obs, 3)
            self.assertEqual(adata.n_vars, 8)
            self.assertIn("mapping_weights", adata.obsm)
            self.assertIn("CXCL13", adata.var_names)
            self.assertIn("MSLN", adata.var_names)
            self.assertEqual(list(adata.obs["top_cell_type"]), ["epithelial", "epithelial", "immune"])


if __name__ == "__main__":
    unittest.main()
