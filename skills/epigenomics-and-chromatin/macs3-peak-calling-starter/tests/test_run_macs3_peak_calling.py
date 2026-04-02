from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"
SCRIPT = ROOT / "skills" / "epigenomics-and-chromatin" / "macs3-peak-calling-starter" / "scripts" / "run_macs3_peak_calling.py"


class Macs3PeakCallingTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_peak(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "macs3_peak_summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--treatment",
                    "skills/epigenomics-and-chromatin/macs3-peak-calling-starter/examples/toy_treatment.bed",
                    "--workdir",
                    str(Path(tmp_dir) / "work"),
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
            self.assertEqual(payload["peak_count"], 1)
            self.assertEqual(payload["treatment_tags"], 20)
            self.assertEqual(payload["fragment_length"], 75)
            self.assertEqual(payload["top_peak"]["name"], "toy_peak_1")
            self.assertEqual(payload["top_peak"]["start"], 102)
            self.assertEqual(payload["top_peak"]["log10_qvalue"], 14.3029)


if __name__ == "__main__":
    unittest.main()
