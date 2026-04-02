from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
AGENTS_PYTHON = ROOT / "slurm" / "envs" / "agents" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-agents-and-automation"
    / "langgraph-planning-execution-agent-starter"
    / "scripts"
    / "run_langgraph_planning_agent.py"
)


def run_goal(goal: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_path = Path(tmp_dir) / "agent.json"
        subprocess.run(
            [
                str(AGENTS_PYTHON),
                str(SCRIPT),
                "--goal",
                goal,
                "--out",
                str(out_path),
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        return json.loads(out_path.read_text(encoding="utf-8"))


class LangGraphPlanningExecutionAgentTests(unittest.TestCase):
    def test_run_agent_routes_marker_and_report_steps(self) -> None:
        payload = run_goal("single-cell marker ranking with an interactive report")
        self.assertEqual(payload["plan_step_count"], 2)
        self.assertIn("scanpy-ranked-genes-starter", payload["unique_skill_slugs"])
        self.assertIn("plotly-interactive-report-starter", payload["unique_skill_slugs"])

    def test_cli_runtime_for_bayesian_publication_goal(self) -> None:
        payload = run_goal("bayesian regression with a publication figure")
        self.assertIn("pymc-bayesian-linear-regression-starter", payload["unique_skill_slugs"])
        self.assertIn("matplotlib-publication-plot-starter", payload["unique_skill_slugs"])


if __name__ == "__main__":
    unittest.main()
