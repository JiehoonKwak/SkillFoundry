from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "genomics" / "fastqc-multiqc-read-qc-starter" / "scripts" / "run_fastqc_multiqc_read_qc.py"


class FastqcMultiqcReadQcTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "read_qc_summary.json"
            workdir = Path(tmp_dir) / "work"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    "skills/genomics/fastqc-multiqc-read-qc-starter/examples/toy_reads.fastq",
                    "--workdir",
                    str(workdir),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_name"], "toy_reads")
            self.assertEqual(payload["total_sequences"], 4)
            self.assertEqual(payload["gc_percent"], 50)
            self.assertEqual(payload["multiqc_sample_count"], 1)
            self.assertTrue(payload["reports"]["multiqc_html"].endswith("multiqc_report.html"))


if __name__ == "__main__":
    unittest.main()
