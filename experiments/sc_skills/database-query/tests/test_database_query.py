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
RUNNER_SCRIPT = SKILL_DIR / "scripts" / "run_database_query.py"
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

            query_rows = read_tsv(outdir / "query_results.tsv")
            normalization_rows = read_tsv(outdir / "normalization_qc.tsv")
            fanout_rows = read_tsv(outdir / "fanout_qc.tsv")
            reconciliation_rows = read_tsv(outdir / "reconciliation_qc.tsv")
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))

            self.assertEqual(len(query_rows), expected["query_count"])
            self.assertEqual(len(normalization_rows), expected["query_count"])
            self.assertEqual(len(reconciliation_rows), expected["query_count"])
            self.assertEqual(len(fanout_rows), 14)
            self.assertEqual(summary["resolved_count"], expected["resolved_count"])
            self.assertEqual(summary["ambiguous_query_count"], len(expected["ambiguous_queries"]))
            self.assertEqual(
                summary["cross_source_supported_query_count"],
                expected["minimum_cross_source_supported_queries"],
            )

            query_index = {row["query"]: row for row in query_rows}
            normalization_index = {row["query"]: row for row in normalization_rows}
            reconciliation_index = {row["query"]: row for row in reconciliation_rows}

            for item in expected["expected_results"]:
                query = item["query"]
                self.assertEqual(query_index[query]["normalized_id"], item["normalized_id"])
                self.assertEqual(
                    normalization_index[query]["inferred_entity_type"],
                    item["inferred_entity_type"],
                )
                observed_sources = set(
                    filter(None, reconciliation_index[query]["supporting_sources"].split(";"))
                )
                self.assertTrue(set(item["supporting_sources"]).issubset(observed_sources))
                self.assertGreater(float(reconciliation_index[query]["score_margin"]), 0.0)

            p53_scores = normalization_index[" p53 "]
            self.assertGreater(
                float(p53_scores["gene_type_score"]),
                float(p53_scores["protein_type_score"]),
            )
            protein_scores = normalization_index["P04637"]
            self.assertGreater(
                float(protein_scores["protein_type_score"]),
                float(protein_scores["gene_type_score"]),
            )
            variant_scores = normalization_index["rs121913343"]
            self.assertGreater(
                float(variant_scores["variant_type_score"]),
                float(variant_scores["gene_type_score"]),
            )
            disease_scores = normalization_index["Li Fraumeni syndrome"]
            self.assertGreater(
                float(disease_scores["disease_type_score"]),
                float(disease_scores["protein_type_score"]),
            )

    def test_runner_fails_for_unresolved_query(self) -> None:
        payload = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
        payload["queries"] = payload["queries"] + ["unknown demo lookup"]

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
            self.assertIn("No hits were found for query", completed.stderr)


if __name__ == "__main__":
    unittest.main()
