from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"
SCRIPT = ROOT / "skills" / "genomics" / "pysam-sam-bam-summary-starter" / "scripts" / "run_pysam_sam_bam_summary.py"
EXAMPLE = ROOT / "skills" / "genomics" / "pysam-sam-bam-summary-starter" / "examples" / "toy_reads.sam"


class PysamSamBamSummaryTests(unittest.TestCase):
    def test_summary_and_bam_write_match_expected_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            bam_path = Path(tmp_dir) / "toy_reads.bam"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path), "--bam-out", str(bam_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["record_count"], 4)
            self.assertEqual(payload["mapped_count"], 3)
            self.assertEqual(payload["unmapped_count"], 1)
            self.assertEqual(payload["paired_count"], 2)
            self.assertEqual(payload["proper_pair_count"], 2)
            self.assertEqual(payload["reference_mapped_counts"], {"chr1": 2, "chr2": 1})
            self.assertAlmostEqual(payload["mean_query_length"], 8.0)
            self.assertTrue(bam_path.exists())
            self.assertTrue(Path(payload["bam_index_path"]).exists())

    def test_rejects_unsupported_extension(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bogus = Path(tmp_dir) / "toy_reads.txt"
            bogus.write_text("not an alignment file\n", encoding="utf-8")
            completed = subprocess.run(
                [str(PYTHON), str(SCRIPT), "--input", str(bogus)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Unsupported input extension", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
