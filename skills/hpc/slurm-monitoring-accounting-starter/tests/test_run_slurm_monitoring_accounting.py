from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "hpc"
    / "slurm-monitoring-accounting-starter"
    / "scripts"
    / "run_slurm_monitoring_accounting.py"
)


class SlurmMonitoringAccountingTests(unittest.TestCase):
    def test_runtime_probe_collects_queue_and_accounting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "slurm_monitoring.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--partition",
                    "cpu",
                    "--job-name",
                    "test-slurm-monitor",
                    "--sleep",
                    "2",
                    "--max-polls",
                    "12",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertGreaterEqual(payload["queue_snapshot_count"], 1)
            self.assertTrue(set(payload["queue_states_seen"]) & {"PENDING", "RUNNING"})
            self.assertEqual(payload["accounting"]["State"], "COMPLETED")
            self.assertEqual(payload["accounting"]["ExitCode"], "0:0")
            self.assertIn("started=", payload["stdout"])
            self.assertEqual(payload["stderr"], "")

    def test_rejects_invalid_poll_interval(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--poll-interval",
                "0",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("poll-interval must be positive", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
