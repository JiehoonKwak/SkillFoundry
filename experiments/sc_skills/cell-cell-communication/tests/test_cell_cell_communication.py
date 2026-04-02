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
    def test_toy_run_writes_expected_outputs_and_numeric_invariants(self) -> None:
        toy_input = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        top_expected = toy_input["expected_invariants"]["top_priority_pair"]

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

            communication_rows = read_tsv(outdir / "communication_results.tsv")
            priority_rows = read_tsv(outdir / "priority_pairs.tsv")
            group_rows = read_tsv(outdir / "qc_group_expression.tsv")
            null_rows = read_tsv(outdir / "qc_pair_nulls.tsv")
            spatial_rows = read_tsv(outdir / "qc_spatial_support.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

            self.assertEqual(priority_rows[0]["ligand"], top_expected["ligand"])
            self.assertEqual(priority_rows[0]["receptor"], top_expected["receptor"])
            self.assertEqual(priority_rows[0]["source_group"], top_expected["source_group"])
            self.assertEqual(priority_rows[0]["target_group"], top_expected["target_group"])
            self.assertGreater(float(priority_rows[0]["score"]), float(priority_rows[1]["score"]))
            self.assertEqual(len(priority_rows), toy_input["priority_top_n"])

            cxcl12_tumor = find_row(group_rows, group="Tumor", gene="CXCL12")
            cxcr4_myeloid = find_row(group_rows, group="Myeloid", gene="CXCR4")
            self.assertAlmostEqual(float(cxcl12_tumor["mean_expression"]), 7.0, places=6)
            self.assertAlmostEqual(float(cxcr4_myeloid["mean_expression"]), 7.333333, places=5)

            null_top = find_row(
                null_rows,
                source_group=top_expected["source_group"],
                target_group=top_expected["target_group"],
                ligand=top_expected["ligand"],
                receptor=top_expected["receptor"],
            )
            self.assertLessEqual(float(null_top["empirical_pvalue"]), 0.125)
            self.assertGreater(float(null_top["observed_mean_stat"]), float(null_top["null_mean_stat"]))

            spatial_top = find_row(
                spatial_rows,
                source_group=top_expected["source_group"],
                target_group=top_expected["target_group"],
                ligand=top_expected["ligand"],
                receptor=top_expected["receptor"],
            )
            self.assertEqual(spatial_top["spatial_support"], "supported")
            self.assertGreaterEqual(
                int(spatial_top["supporting_edges"]),
                toy_input["expected_invariants"]["positive_pair_min_spatial_edges"],
            )
            self.assertAlmostEqual(float(spatial_top["adjacency_fraction"]), 1.0, places=6)

            top_result = communication_rows[0]
            self.assertEqual(top_result["ligand"], top_expected["ligand"])
            self.assertEqual(top_result["receptor"], top_expected["receptor"])
            self.assertEqual(top_result["spatial_support"], "supported")
            self.assertGreater(float(top_result["score"]), 0.8)

            self.assertEqual(summary["top_pair"]["ligand"], top_expected["ligand"])
            self.assertEqual(summary["counts"]["communication_rows"], len(communication_rows))
            self.assertIn("interpretation_report.md", summary["written_files"])


if __name__ == "__main__":
    unittest.main()
