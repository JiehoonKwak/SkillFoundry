from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NUMERICS_PYTHON = ROOT / "slurm" / "envs" / "numerics" / "bin" / "python"
AGENTS_PYTHON = ROOT / "slurm" / "envs" / "agents" / "bin" / "python"
INSTRUMENTATION_PYTHON = ROOT / "slurm" / "envs" / "instrumentation" / "bin" / "python"
REPORTING_PYTHON = ROOT / "slurm" / "envs" / "reporting" / "bin" / "python"


class Phase38FrontierConversionSmokeTests(unittest.TestCase):
    def test_chaospy_uncertainty_propagation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(NUMERICS_PYTHON),
                    "skills/scientific-computing-and-numerical-methods/chaospy-uncertainty-propagation-starter/scripts/run_chaospy_uncertainty_propagation.py",
                    "--config",
                    "skills/scientific-computing-and-numerical-methods/chaospy-uncertainty-propagation-starter/examples/toy_parameters.json",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["node_count"], 16)
            self.assertAlmostEqual(payload["mean"], 1.51, places=6)

    def test_inspect_evaluation_harness_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            out_path = tmp_root / "summary.json"
            log_dir = tmp_root / "logs"
            subprocess.run(
                [
                    str(AGENTS_PYTHON),
                    "skills/scientific-agents-and-automation/inspect-evaluation-harness-starter/scripts/run_inspect_evaluation_harness.py",
                    "--cases",
                    "skills/scientific-agents-and-automation/inspect-evaluation-harness-starter/examples/toy_eval_cases.json",
                    "--summary-out",
                    str(out_path),
                    "--log-dir",
                    str(log_dir),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertAlmostEqual(payload["candidate_accuracy"], 1.0, places=6)
            self.assertAlmostEqual(payload["baseline_accuracy"], 0.25, places=6)

    def test_autoprotocol_experiment_plan_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            protocol_path = tmp_root / "protocol.json"
            summary_path = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(INSTRUMENTATION_PYTHON),
                    "skills/robotics-lab-automation-and-scientific-instrumentation/autoprotocol-experiment-plan-starter/scripts/build_autoprotocol_experiment_plan.py",
                    "--transfers",
                    "skills/robotics-lab-automation-and-scientific-instrumentation/autoprotocol-experiment-plan-starter/examples/toy_transfers.tsv",
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
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["transfer_count"], 2)
            self.assertIn("incubate", payload["instruction_ops"])

    def test_dash_scientific_dashboard_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            html_out = tmp_root / "dashboard_preview.html"
            summary_out = tmp_root / "dashboard_summary.json"
            subprocess.run(
                [
                    str(REPORTING_PYTHON),
                    "skills/visualization-and-reporting/dash-scientific-dashboard-starter/scripts/build_dash_scientific_dashboard.py",
                    "--input",
                    "skills/visualization-and-reporting/dash-scientific-dashboard-starter/examples/toy_measurements.tsv",
                    "--html-out",
                    str(html_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["callback_count"], 1)
            self.assertEqual(payload["trace_count"], 2)


if __name__ == "__main__":
    unittest.main()
