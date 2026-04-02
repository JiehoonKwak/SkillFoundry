from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "instrumentation" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "robotics-lab-automation-and-scientific-instrumentation"
    / "autoprotocol-experiment-plan-starter"
    / "scripts"
    / "build_autoprotocol_experiment_plan.py"
)
EXAMPLES = (
    ROOT
    / "skills"
    / "robotics-lab-automation-and-scientific-instrumentation"
    / "autoprotocol-experiment-plan-starter"
    / "examples"
)


class AutoprotocolExperimentPlanTests(unittest.TestCase):
    def test_protocol_and_summary_are_generated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            protocol_path = tmp_root / "protocol.json"
            summary_path = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--transfers",
                    str(EXAMPLES / "toy_transfers.tsv"),
                    "--protocol-json",
                    str(protocol_path),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            protocol_payload = json.loads(protocol_path.read_text(encoding="utf-8"))
            summary_payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(sorted(protocol_payload["refs"].keys()), ["assay_plate", "source_plate"])
            self.assertEqual(summary_payload["transfer_count"], 2)
            self.assertIn("liquid_handle", summary_payload["instruction_ops"])
            self.assertIn("incubate", summary_payload["instruction_ops"])
            self.assertTrue(summary_payload["sealed_before_incubate"])

    def test_missing_transfer_plan_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--transfers",
                "missing.tsv",
                "--protocol-json",
                "ignored.json",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Transfer plan not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
