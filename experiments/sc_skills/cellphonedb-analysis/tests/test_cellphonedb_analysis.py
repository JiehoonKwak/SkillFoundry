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
    def test_toy_run_computes_cpdb_like_statistics_and_qc(self) -> None:
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

            significant_rows = read_tsv(outdir / "cellphonedb_significant_means.tsv")
            ranked_rows = read_tsv(outdir / "interaction_ranked.tsv")
            group_rows = read_tsv(outdir / "qc_group_expression.tsv")
            null_rows = read_tsv(outdir / "qc_statistical_null.tsv")
            deg_rows = read_tsv(outdir / "qc_deg_support.tsv")
            adjacency_rows = read_tsv(outdir / "qc_adjacency_support.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

            self.assertGreaterEqual(
                len(significant_rows),
                toy_input["expected_invariants"]["minimum_significant_rows"],
            )
            self.assertEqual(significant_rows[0]["interacting_pair"], expected_top["interacting_pair"])
            self.assertEqual(significant_rows[0]["source_group"], expected_top["source_group"])
            self.assertEqual(significant_rows[0]["target_group"], expected_top["target_group"])
            self.assertLessEqual(float(significant_rows[0]["pvalue"]), toy_input["pvalue_cutoff"])

            cxcl12_tumor = find_row(group_rows, group="Tumor", gene="CXCL12")
            cxcr4_myeloid = find_row(group_rows, group="Myeloid", gene="CXCR4")
            self.assertAlmostEqual(float(cxcl12_tumor["mean_expression"]), 7.0, places=6)
            self.assertAlmostEqual(float(cxcr4_myeloid["mean_expression"]), 7.333333, places=5)

            top_null = find_row(
                null_rows,
                interacting_pair=expected_top["interacting_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertLessEqual(float(top_null["empirical_pvalue"]), toy_input["pvalue_cutoff"])
            self.assertGreater(
                float(top_null["observed_mean_expression"]),
                float(top_null["null_mean_expression"]),
            )

            top_deg = find_row(
                deg_rows,
                interacting_pair=expected_top["interacting_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_deg["ligand_is_deg"], "true")
            self.assertEqual(top_deg["receptor_is_deg"], "true")
            self.assertEqual(top_deg["deg_relevant"], "true")

            top_adjacency = find_row(
                adjacency_rows,
                interacting_pair=expected_top["interacting_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_adjacency["spatial_support"], "supported")
            self.assertGreaterEqual(
                int(top_adjacency["supporting_edges"]),
                toy_input["expected_invariants"]["positive_pair_min_spatial_edges"],
            )
            self.assertAlmostEqual(float(top_adjacency["adjacency_fraction"]), 1.0, places=6)

            self.assertEqual(ranked_rows[0]["interacting_pair"], expected_top["interacting_pair"])
            self.assertEqual(ranked_rows[0]["supporting_group_pair"], "Tumor->Myeloid")
            self.assertIn("significant_mean", ranked_rows[0]["evidence"])
            self.assertGreater(float(ranked_rows[0]["score"]), float(ranked_rows[1]["score"]))

            self.assertEqual(summary["top_interaction"]["interacting_pair"], expected_top["interacting_pair"])
            self.assertEqual(summary["counts"]["significant_rows"], len(significant_rows))
            self.assertIn("cellphonedb_report.md", summary["written_files"])


if __name__ == "__main__":
    unittest.main()
