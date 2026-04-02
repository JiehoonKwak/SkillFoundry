from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_latent_qc_and_label_transfer(self) -> None:
        payload = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        expected = payload["expected_invariants"]

        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            completed = subprocess.run(
                ["python3", str(SCRIPT), "--outdir", str(outdir)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            qc = pd.read_csv(outdir / "modality_qc.tsv", sep="\t")
            qc_lookup = {
                (row.metric, row.modality): float(row.value)
                for row in qc.itertuples(index=False)
            }
            self.assertGreaterEqual(qc_lookup[("batch_mixing", "joint")], expected["min_batch_mixing"])
            self.assertGreaterEqual(qc_lookup[("label_consistency", "joint")], expected["min_label_consistency"])
            self.assertAlmostEqual(qc_lookup[("query_label_accuracy", "joint")], expected["min_query_accuracy"], places=6)

            for modality in ("rna", "protein", "chromatin"):
                self.assertGreater(
                    qc_lookup[("between_label_centroid_distance", modality)],
                    qc_lookup[("between_batch_centroid_distance", modality)],
                )

            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["counts"]["cells"], 8)
            predictions = {row["cell_id"]: row for row in summary["query_predictions"]}
            self.assertEqual(
                {cell_id: row["predicted_label"] for cell_id, row in predictions.items()},
                expected["query_labels"],
            )
            self.assertLess(
                predictions[expected["ambiguous_query_cell"]]["confidence"],
                predictions["cell_005"]["confidence"],
            )

            report_text = (outdir / "integration_report.md").read_text(encoding="utf-8")
            self.assertIn("Model choice: multivi_style", report_text)
            self.assertIn("## Label transfer", report_text)

            with h5py.File(outdir / "integrated_latent.h5mu", "r") as handle:
                self.assertEqual(handle.attrs["format"], "mudata_like_h5mu_surrogate")
                self.assertEqual(handle["mod"]["rna"]["raw_counts"].shape, (8, 4))
                self.assertEqual(handle["mod"]["protein"]["raw_counts"].shape, (8, 2))
                self.assertEqual(handle["mod"]["chromatin"]["raw_counts"].shape, (8, 3))
                self.assertEqual(handle["obsm"]["X_latent"].shape, (8, 2))
                self.assertEqual(handle["obsp"]["knn_indices"].shape, (8, 3))
                query_ids = [
                    value.decode("utf-8") if isinstance(value, bytes) else value
                    for value in handle["label_transfer"]["query_ids"][()]
                ]
                predicted_labels = [
                    value.decode("utf-8") if isinstance(value, bytes) else value
                    for value in handle["label_transfer"]["predicted_label"][()]
                ]
                self.assertEqual(dict(zip(query_ids, predicted_labels, strict=True)), expected["query_labels"])


if __name__ == "__main__":
    unittest.main()
