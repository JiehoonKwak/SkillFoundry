from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
EXERCISE_SCRIPT = SKILL_DIR / "scripts" / "run_exercise.py"
RUNNER_SCRIPT = SKILL_DIR / "scripts" / "run_panel_design.py"
INPUT_PATH = SKILL_DIR / "examples" / "toy_input.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SkillLocalToyRunTests(unittest.TestCase):
    def test_toy_run_computes_expected_invariants(self) -> None:
        payload = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
        expected = payload["expected_invariants"]

        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            completed = subprocess.run(
                ["python3", str(EXERCISE_SCRIPT), "--outdir", str(outdir)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            panel_rows = read_tsv(outdir / "panel_candidates.tsv")
            aggregation_rows = read_tsv(outdir / "aggregated_marker_qc.tsv")
            redundancy_rows = read_tsv(outdir / "redundancy_qc.tsv")
            coverage_rows = read_tsv(outdir / "coverage_balance_qc.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

            self.assertEqual(len(aggregation_rows), expected["aggregated_candidate_count"])
            self.assertEqual(summary["selected_gene_count"], expected["selected_gene_count"])
            self.assertEqual(len(panel_rows), expected["selected_gene_count"])
            self.assertEqual(len(coverage_rows), expected["selection_trace_steps"])
            self.assertEqual(len({row["panel_role"] for row in panel_rows}), expected["minimum_role_count"])
            self.assertGreaterEqual(summary["retained_candidate_count"], expected["minimum_retained_candidates"])

            selected_gene_set = {row["target_gene"] for row in panel_rows}
            self.assertEqual(selected_gene_set, set(expected["selected_gene_set"]))

            pruned_gene_set = {
                row["gene"] for row in redundancy_rows if row["retained_after_pruning"] == "no"
            }
            self.assertEqual(pruned_gene_set, set(expected["pruned_gene_set"]))
            self.assertEqual(
                len({row["redundancy_component"] for row in redundancy_rows}),
                expected["redundancy_component_count"],
            )

            ranking = {int(row["rank"]): row["target_gene"] for row in panel_rows}
            self.assertEqual(ranking[1], expected["top_ranked_gene"])
            self.assertIn(expected["flex_gene"], [row["selected_gene"] for row in coverage_rows])

            first_step = coverage_rows[0]
            self.assertEqual(first_step["selection_stage"], "quota")
            self.assertEqual(first_step["selected_gene"], "EPCAM")

            final_step = coverage_rows[-1]
            self.assertEqual(final_step["selection_stage"], "balance_fill")
            self.assertEqual(final_step["selected_gene"], expected["flex_gene"])
            self.assertGreater(float(final_step["marginal_gain"]), 0.1)

            for cell_type, minimum_coverage in expected["minimum_final_coverage"].items():
                self.assertGreaterEqual(summary["coverage_by_cell_type"][cell_type], minimum_coverage)

            redundancy_index = {row["gene"]: row for row in redundancy_rows}
            self.assertEqual(redundancy_index["KRT19"]["representative_gene"], "EPCAM")
            self.assertEqual(redundancy_index["TRBC1"]["representative_gene"], "CD3D")
            self.assertEqual(redundancy_index["VWF"]["representative_gene"], "KDR")

    def test_runner_fails_for_infeasible_panel_size(self) -> None:
        payload = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
        payload["panel_size"] = 4

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            bad_input = tmpdir_path / "bad_input.json"
            outdir = tmpdir_path / "out"
            bad_input.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

            completed = subprocess.run(
                ["python3", str(RUNNER_SCRIPT), "--input", str(bad_input), "--outdir", str(outdir)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Panel size is smaller than the required cell-type quotas.", completed.stderr)


if __name__ == "__main__":
    unittest.main()
