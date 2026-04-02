from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "metadata-harmonization-starter"
SCRIPT = SKILL / "scripts" / "run_metadata_harmonization.py"
COHORT_A = SKILL / "examples" / "cohort_a.tsv"
COHORT_B = SKILL / "examples" / "cohort_b.tsv"
MAPPING = SKILL / "examples" / "column_mapping.json"


class MetadataHarmonizationTests(unittest.TestCase):
    def run_skill(self, *extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *extra],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_harmonizes_two_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_tsv = Path(tmpdir) / "harmonized.tsv"
            summary = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(COHORT_A),
                "--input", str(COHORT_B),
                "--mapping", str(MAPPING),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(summary),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["row_count"], 4)
            self.assertEqual(payload["source_counts"], {"cohort_a.tsv": 2, "cohort_b.tsv": 2})
            self.assertEqual(payload["distinct_conditions"], ["control", "treated"])
            self.assertEqual(payload["first_row"]["sample_id"], "A1")

    def test_rejects_unmapped_input_basename(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            unknown = Path(tmpdir) / "unknown.tsv"
            unknown.write_text("id\tvalue\nX\t1\n", encoding="utf-8")
            out_tsv = Path(tmpdir) / "harmonized.tsv"
            summary = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(unknown),
                "--mapping", str(MAPPING),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(summary),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No source mapping", result.stderr)

    def test_rejects_missing_mapped_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            broken = Path(tmpdir) / "cohort_a.tsv"
            broken.write_text("sample_id\tdonor\tsex\tage_years\nA1\tD100\tF\t34\n", encoding="utf-8")
            out_tsv = Path(tmpdir) / "harmonized.tsv"
            summary = Path(tmpdir) / "summary.json"
            result = self.run_skill(
                "--input", str(broken),
                "--mapping", str(MAPPING),
                "--out-tsv", str(out_tsv),
                "--summary-out", str(summary),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("missing mapped field", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
