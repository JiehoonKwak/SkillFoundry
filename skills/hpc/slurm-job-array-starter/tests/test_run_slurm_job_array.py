from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = ROOT / "skills" / "hpc" / "slurm-job-array-starter" / "scripts" / "run_slurm_job_array.py"
SPEC = importlib.util.spec_from_file_location("slurm_job_array", SCRIPT_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SlurmJobArrayTests(unittest.TestCase):
    def test_expand_array_spec_supports_ranges_and_steps(self) -> None:
        self.assertEqual(MODULE.expand_array_spec("0-4:2,7"), [0, 2, 4, 7])

    def test_runtime_summary_reports_completed_tasks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "array_summary.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT_PATH),
                    "--partition",
                    "cpu",
                    "--job-name",
                    "slurm-array-test",
                    "--array-spec",
                    "0-1",
                    "--sleep",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["task_count"], 2)
            self.assertEqual(payload["task_ids"], [0, 1])
            self.assertGreaterEqual(payload["queue_snapshot_count"], 1)
            for task in payload["task_summaries"]:
                self.assertEqual(task["accounting"]["State"], "COMPLETED")
                self.assertEqual(task["accounting"]["ExitCode"], "0:0")
                self.assertIn(f"array_task_id={task['task_id']}", task["stdout"])
                self.assertEqual(task["stderr"], "")


if __name__ == "__main__":
    unittest.main()
