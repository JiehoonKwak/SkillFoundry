from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MINIWDL = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "miniwdl"
UDOCKER = ROOT / "slurm" / "envs" / "workflow-languages" / "bin" / "udocker"


class WdlTaskStarterTests(unittest.TestCase):
    @unittest.skipUnless(MINIWDL.exists() and UDOCKER.exists(), f"miniwdl or udocker missing under {ROOT / 'slurm' / 'envs' / 'workflow-languages' / 'bin'}")
    def test_run_wdl_hello_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/wdl-task-starter/scripts/run_wdl_hello.py",
                    "--name",
                    "WDL unittest",
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
            self.assertEqual(payload["output_text"], "Hello, WDL unittest")
            self.assertEqual(payload["greeting_text"], "Hello, WDL unittest")
            self.assertIn("hello_workflow.greeting_file", payload["outputs"]["outputs"])


if __name__ == "__main__":
    unittest.main()
