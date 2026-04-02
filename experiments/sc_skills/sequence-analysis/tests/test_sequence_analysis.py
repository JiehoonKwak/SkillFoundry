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
RUNNER_SCRIPT = SKILL_DIR / "scripts" / "run_sequence_analysis.py"
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

            sequence_rows = read_tsv(outdir / "sequence_summary.tsv")
            feature_rows = read_tsv(outdir / "feature_annotations.tsv")
            primer_rows = read_tsv(outdir / "primer_candidates.tsv")
            lookup_rows = read_tsv(outdir / "lookup_qc.tsv")
            variant_rows = read_tsv(outdir / "variant_qc.tsv")
            primer_qc_rows = read_tsv(outdir / "primer_qc.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

            self.assertEqual(len(sequence_rows), expected["query_count"])
            self.assertEqual(len(feature_rows), expected["minimum_feature_annotation_count"])
            self.assertEqual(len(primer_rows), expected["query_count"] * expected["primer_pairs_per_query"])
            self.assertEqual(len(lookup_rows), expected["query_count"] * expected["lookup_candidates_per_query"])
            self.assertEqual(len(variant_rows), expected["query_count"])
            self.assertGreater(len(primer_qc_rows), len(primer_rows))

            self.assertEqual(summary["query_count"], expected["query_count"])
            self.assertEqual(summary["resolved_query_count"], expected["query_count"])
            self.assertEqual(
                summary["reverse_complement_query_count"],
                expected["reverse_complement_query_count"],
            )
            self.assertEqual(summary["mutant_variant_count"], expected["mutant_variant_count"])
            self.assertEqual(summary["feature_annotation_count"], expected["minimum_feature_annotation_count"])
            self.assertEqual(summary["primer_pair_count"], len(primer_rows))
            self.assertIn("run_summary.json", summary["written_files"])

            sequence_index = {row["query_id"]: row for row in sequence_rows}
            variant_index = {row["query_id"]: row for row in variant_rows}
            top_lookup = {row["query_id"]: row for row in lookup_rows if row["candidate_rank"] == "1"}
            top_primer = {row["query_id"]: row for row in primer_rows if row["pair_rank"] == "1"}

            for query_id, matched_symbol in expected["matched_symbols"].items():
                self.assertEqual(sequence_index[query_id]["matched_symbol"], matched_symbol)
                self.assertEqual(sequence_index[query_id]["matched_accession"], expected["matched_accessions"][query_id])
                self.assertEqual(sequence_index[query_id]["orientation"], expected["orientation_by_query"][query_id])
                self.assertGreaterEqual(
                    float(sequence_index[query_id]["identity_pct"]),
                    expected["minimum_identity_pct"],
                )
                self.assertEqual(variant_index[query_id]["observed_state"], expected["variant_state_by_query"][query_id])
                self.assertEqual(top_lookup[query_id]["matched_symbol"], matched_symbol)
                self.assertEqual(
                    int(top_primer[query_id]["product_size"]),
                    expected["top_product_size_by_query"][query_id],
                )
                self.assertEqual(
                    top_primer[query_id]["amplicon_target_label"],
                    expected["amplicon_target_by_query"][query_id],
                )

            self.assertEqual(top_lookup["query_001"]["orientation"], "+")
            self.assertEqual(top_lookup["query_002"]["orientation"], "-")
            self.assertGreater(
                float(top_lookup["query_001"]["score"]),
                float([row for row in lookup_rows if row["query_id"] == "query_001" and row["candidate_rank"] == "2"][0]["score"]),
            )
            self.assertGreater(
                float(top_lookup["query_002"]["score"]),
                float([row for row in lookup_rows if row["query_id"] == "query_002" and row["candidate_rank"] == "2"][0]["score"]),
            )

            hotspot_rows = {
                (row["query_id"], row["feature_id"])
                for row in feature_rows
                if row["feature_type"] == "hotspot_region"
            }
            self.assertEqual(
                hotspot_rows,
                {("query_001", "BRAF_hotspot_window"), ("query_002", "KRAS_hotspot_window")},
            )

            product_min, product_max = payload["primer_constraints"]["product_size_range"]
            for row in primer_rows:
                product_size = int(row["product_size"])
                self.assertGreaterEqual(product_size, product_min)
                self.assertLessEqual(product_size, product_max)
                self.assertLessEqual(abs(float(row["left_tm"]) - float(row["right_tm"])), payload["primer_constraints"]["max_tm_gap"])
                self.assertGreaterEqual(float(row["left_gc_pct"]), payload["primer_constraints"]["gc_range"][0])
                self.assertLessEqual(float(row["right_gc_pct"]), payload["primer_constraints"]["gc_range"][1])

    def test_runner_fails_when_primer_window_is_infeasible(self) -> None:
        payload = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
        payload["primer_constraints"]["product_size_range"] = [140, 150]

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
            self.assertIn("No acceptable primer pairs found", completed.stderr)


if __name__ == "__main__":
    unittest.main()
