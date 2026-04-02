from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_row(rows: list[dict[str, str]], **filters: str) -> dict[str, str]:
    for row in rows:
        if all(row[key] == value for key, value in filters.items()):
            return row
    raise AssertionError(f"Missing row matching filters: {filters}")


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_liana_like_ranking_and_qc(self) -> None:
        toy_input = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        expected_top = toy_input["expected_invariants"]["top_interaction"]

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

            rankings = read_tsv(outdir / "liana_rankings.tsv")
            resource_rows = read_tsv(outdir / "resource_overlap.tsv")
            harmonization_rows = read_tsv(outdir / "qc_metadata_harmonization.tsv")
            group_rows = read_tsv(outdir / "qc_group_expression.tsv")
            method_rows = read_tsv(outdir / "qc_method_scores.tsv")
            condition_rows = read_tsv(outdir / "qc_condition_summary.tsv")
            adjacency_rows = read_tsv(outdir / "qc_adjacency_support.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            report_text = (outdir / "liana_report.md").read_text(encoding="utf-8")

            self.assertGreaterEqual(len(rankings), toy_input["expected_invariants"]["minimum_ranked_rows"])
            self.assertEqual(rankings[0]["pair"], expected_top["pair"])
            self.assertEqual(rankings[0]["source_group"], expected_top["source_group"])
            self.assertEqual(rankings[0]["target_group"], expected_top["target_group"])
            self.assertGreater(float(rankings[0]["aggregate_score"]), float(rankings[1]["aggregate_score"]))
            self.assertEqual(rankings[0]["dominant_condition"], "treated")
            self.assertEqual(rankings[0]["selected_methods"], ",".join(toy_input["selected_methods"]))

            top_resource = find_row(resource_rows, pair=expected_top["pair"])
            self.assertEqual(
                int(top_resource["resource_count"]),
                toy_input["expected_invariants"]["positive_pair_resource_count"],
            )
            self.assertIn("CellPhoneDB", top_resource["resources"])
            self.assertIn("OmniPath", top_resource["resources"])
            self.assertIn("LIANA", top_resource["resources"])

            self.assertEqual(len({row["normalized_group"] for row in harmonization_rows}), toy_input["expected_invariants"]["group_count"])
            self.assertIn("true", {row["harmonized"] for row in harmonization_rows})

            tumor_cxcl12 = find_row(group_rows, group="Tumor", condition="all", gene="CXCL12")
            myeloid_cxcr4 = find_row(group_rows, group="Myeloid", condition="all", gene="CXCR4")
            tumor_treated_cxcl12 = find_row(group_rows, group="Tumor", condition="treated", gene="CXCL12")
            self.assertAlmostEqual(float(tumor_cxcl12["mean_expression"]), 7.0, places=6)
            self.assertAlmostEqual(float(myeloid_cxcr4["mean_expression"]), 7.0, places=6)
            self.assertAlmostEqual(float(tumor_treated_cxcl12["mean_expression"]), 7.5, places=6)

            top_method = find_row(
                method_rows,
                pair=expected_top["pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_method["expression_pass"], "true")
            self.assertLessEqual(float(top_method["empirical_pvalue"]), 0.1)
            self.assertGreater(float(top_method["geometric_mean"]), float(top_method["null_mean_score"]))
            self.assertAlmostEqual(float(top_method["min_expr"]), 7.0, places=6)

            top_condition = find_row(
                condition_rows,
                pair=toy_input["expected_invariants"]["top_condition_gain_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_condition["dominant_condition"], "treated")
            self.assertGreater(float(top_condition["delta_geometric_mean"]), 1.0)
            self.assertAlmostEqual(float(top_condition["baseline_geometric_mean"]), 6.0, places=6)
            self.assertAlmostEqual(float(top_condition["comparison_geometric_mean"]), 7.5, places=6)

            top_adjacency = find_row(
                adjacency_rows,
                pair=expected_top["pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_adjacency["spatial_support"], "supported")
            self.assertGreaterEqual(
                int(top_adjacency["supporting_edges"]),
                toy_input["expected_invariants"]["positive_pair_min_supporting_edges"],
            )
            self.assertGreater(float(top_adjacency["adjacency_fraction"]), 0.3)

            self.assertEqual(summary["top_interaction"]["pair"], expected_top["pair"])
            self.assertEqual(summary["top_condition_delta"]["pair"], toy_input["expected_invariants"]["top_condition_gain_pair"])
            self.assertEqual(summary["counts"]["cells"], toy_input["expected_invariants"]["cell_count"])
            self.assertEqual(summary["selected_methods"], toy_input["selected_methods"])
            self.assertIn("qc_method_scores.tsv", summary["written_files"])
            self.assertIn("## Condition-aware summary", report_text)


if __name__ == "__main__":
    unittest.main()
