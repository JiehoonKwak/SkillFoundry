from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CWLTOOL = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "cwltool"


class CwlCommandLineToolStarterTests(unittest.TestCase):
    @unittest.skipUnless(CWLTOOL.exists(), f"cwltool not available at {CWLTOOL}")
    def test_run_cwl_hello_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/cwl-commandlinetool-starter/scripts/run_cwl_hello.py",
                    "--message",
                    "hello from unittest",
                    "--workspace",
                    str(Path(tmp_dir) / "workspace"),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["output_text"], "hello from unittest")
            self.assertEqual(payload["raw_summary"]["greeting_file"]["basename"], "greeting.txt")


if __name__ == "__main__":
    unittest.main()
