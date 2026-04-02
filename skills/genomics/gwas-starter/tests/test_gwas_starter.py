from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "genomics" / "gwas-starter"
SCRIPT = SKILL / "scripts" / "run_gwas_summary_qc.py"
TOY = SKILL / "examples" / "toy_sumstats.tsv"
CONFIG = SKILL / "examples" / "qc_config.json"


class GwasSummaryStatisticsStarterTests(unittest.TestCase):
    def run_skill(self, *extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *extra],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )

    def test_qc_summary_and_flagged_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_tsv = Path(tmpdir) / "flagged.tsv"
            out_json = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(TOY),
                "--config", str(CONFIG),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(out_json),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["row_count"], 7)
            self.assertEqual(payload["qc_pass_count"], 3)
            self.assertEqual(payload["qc_fail_count"], 4)
            self.assertEqual(payload["genome_wide_significant_count"], 1)
            self.assertEqual(payload["lead_variants"][0]["variant_id"], "rs111111")
            self.assertEqual(payload["flag_counts"]["duplicate_variant_id"], 1)
            self.assertEqual(payload["recommended_resources"][0]["resource_id"], "gwaslab-docs")

            with out_tsv.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 7)
            synthesized = [row for row in rows if row["variant_id"] == "1:10200:C:T"]
            self.assertEqual(len(synthesized), 1)
            self.assertEqual(synthesized[0]["qc_status"], "pass")
            low_info = [row for row in rows if row["variant_id"] == "rsLOWINFO"][0]
            self.assertIn("low_info", low_info["qc_flags"])

    def test_accepts_odds_ratio_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "odds_ratio_sumstats.tsv"
            input_path.write_text(
                "\n".join(
                    [
                        "variant\tchrom\tposition\teffect_allele\tother_allele\tOR\tSE\tP\tN\tEAF\tINFO",
                        "rsALT\t7\t70101\tA\tG\t1.25\t0.05\t2e-7\t42000\t0.22\t0.97",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_tsv = Path(tmpdir) / "flagged.tsv"
            out_json = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(input_path),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(out_json),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["qc_pass_count"], 1)
            self.assertEqual(payload["lead_variants"][0]["variant_id"], "rsALT")
            with out_tsv.open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertNotEqual(row["beta"], "")

    def test_rejects_missing_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "broken.tsv"
            input_path.write_text(
                "SNP\tCHR\tEA\tNEA\tBETA\tSE\tP\nrs1\t1\tA\tG\t0.1\t0.01\t0.05\n",
                encoding="utf-8",
            )
            out_tsv = Path(tmpdir) / "flagged.tsv"
            out_json = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(input_path),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(out_json),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Missing required GWAS columns", result.stderr)


if __name__ == "__main__":
    unittest.main()
