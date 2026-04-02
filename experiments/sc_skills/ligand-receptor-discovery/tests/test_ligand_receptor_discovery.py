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
    def test_toy_run_computes_ranked_candidates_and_qc(self) -> None:
        toy_input = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        expected_top = toy_input["expected_invariants"]["top_candidate"]

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

            candidate_rows = read_tsv(outdir / "candidate_pairs.tsv")
            evidence_rows = read_tsv(outdir / "evidence_table.tsv")
            group_rows = read_tsv(outdir / "qc_group_expression.tsv")
            pair_rows = read_tsv(outdir / "qc_pair_scores.tsv")
            null_rows = read_tsv(outdir / "qc_shuffle_null.tsv")
            spatial_rows = read_tsv(outdir / "qc_spatial_support.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            report_text = (outdir / "discovery_summary.md").read_text(encoding="utf-8")

            self.assertGreaterEqual(len(candidate_rows), toy_input["expected_invariants"]["minimum_candidate_rows"])
            self.assertEqual(candidate_rows[0]["candidate_pair"], expected_top["candidate_pair"])
            self.assertEqual(candidate_rows[0]["source_group"], expected_top["source_group"])
            self.assertEqual(candidate_rows[0]["target_group"], expected_top["target_group"])
            self.assertGreater(float(candidate_rows[0]["support_score"]), float(candidate_rows[1]["support_score"]))
            self.assertEqual(candidate_rows[0]["spatial_support"], "supported")

            tumor_cxcl12 = find_row(group_rows, group="Tumor", gene="CXCL12")
            myeloid_cxcr4 = find_row(group_rows, group="Myeloid", gene="CXCR4")
            self.assertAlmostEqual(float(tumor_cxcl12["mean_expression"]), 7.0, places=6)
            self.assertAlmostEqual(float(tumor_cxcl12["expression_fraction"]), 1.0, places=6)
            self.assertAlmostEqual(float(myeloid_cxcr4["mean_expression"]), 7.0, places=6)

            top_pair_qc = find_row(
                pair_rows,
                candidate_pair=expected_top["candidate_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_pair_qc["expression_pass"], "true")
            self.assertAlmostEqual(float(top_pair_qc["min_expr"]), 7.0, places=6)
            self.assertAlmostEqual(float(top_pair_qc["geometric_mean"]), 7.0, places=6)
            self.assertEqual(
                int(top_pair_qc["resource_count"]),
                toy_input["expected_invariants"]["positive_pair_resource_count"],
            )
            self.assertAlmostEqual(float(top_pair_qc["knowledge_score"]), 1.0, places=6)
            self.assertIn("CellPhoneDB", top_pair_qc["resources"])
            self.assertIn("OmniPath", top_pair_qc["resources"])
            self.assertIn("LIANA", top_pair_qc["resources"])

            top_null = find_row(
                null_rows,
                candidate_pair=expected_top["candidate_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertLessEqual(float(top_null["empirical_pvalue"]), 0.1)
            self.assertGreater(float(top_null["observed_expression_score"]), float(top_null["null_mean_score"]))

            top_spatial = find_row(
                spatial_rows,
                candidate_pair=expected_top["candidate_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
            )
            self.assertEqual(top_spatial["spatial_support"], "supported")
            self.assertGreaterEqual(
                int(top_spatial["supporting_edges"]),
                toy_input["expected_invariants"]["positive_pair_min_supporting_edges"],
            )
            self.assertAlmostEqual(float(top_spatial["adjacency_fraction"]), 1.0, places=6)

            resource_evidence = find_row(
                evidence_rows,
                candidate_pair=expected_top["candidate_pair"],
                source_group=expected_top["source_group"],
                target_group=expected_top["target_group"],
                evidence_type="resources",
            )
            self.assertIn("CellPhoneDB", resource_evidence["evidence_value"])
            self.assertIn("OmniPath", resource_evidence["evidence_value"])

            self.assertEqual(summary["top_candidate"]["candidate_pair"], expected_top["candidate_pair"])
            self.assertEqual(summary["counts"]["cells"], toy_input["expected_invariants"]["cell_count"])
            self.assertIn("candidate_pairs.tsv", summary["written_files"])
            self.assertIn("## Evidence summary", report_text)


if __name__ == "__main__":
    unittest.main()
