from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AGENTS_PYTHON = ROOT / "slurm" / "envs" / "agents" / "bin" / "python"
STATISTICS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class Phase25AgentClinicalProteomicsGraphSmokeTests(unittest.TestCase):
    def test_langgraph_planning_execution_agent_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "langgraph_agent.json"
            subprocess.run(
                [
                    str(AGENTS_PYTHON),
                    "skills/scientific-agents-and-automation/langgraph-planning-execution-agent-starter/scripts/run_langgraph_planning_agent.py",
                    "--goal",
                    "single-cell marker ranking with an interactive report",
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
            self.assertEqual(payload["plan_step_count"], 2)
            self.assertIn("scanpy-ranked-genes-starter", payload["unique_skill_slugs"])

    def test_lifelines_kaplan_meier_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_out = Path(tmp_dir) / "lifelines.json"
            png_out = Path(tmp_dir) / "lifelines.png"
            subprocess.run(
                [
                    str(STATISTICS_PYTHON),
                    "skills/clinical-biomedical-data-science/lifelines-kaplan-meier-starter/scripts/run_lifelines_kaplan_meier.py",
                    "--input",
                    "skills/clinical-biomedical-data-science/lifelines-kaplan-meier-starter/examples/toy_survival_cohort.tsv",
                    "--summary-out",
                    str(summary_out),
                    "--png-out",
                    str(png_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["median_survival_by_group"]["A"], 5.0)
            self.assertTrue(png_out.exists())

    def test_interpro_entry_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "interpro.json"
            subprocess.run(
                [
                    "python3",
                    "skills/proteomics/interpro-entry-summary/scripts/fetch_interpro_entry.py",
                    "--accession",
                    "IPR000023",
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
            self.assertEqual(payload["accession"], "IPR000023")
            self.assertEqual(payload["type"], "domain")

    def test_networkx_graph_construction_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "networkx.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/networkx-graph-construction-starter/scripts/run_networkx_graph_construction.py",
                    "--input",
                    "skills/systems-biology/networkx-graph-construction-starter/examples/toy_pathway_edges.tsv",
                    "--source-node",
                    "EGFR",
                    "--target-node",
                    "STAT3",
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
            self.assertEqual(payload["node_count"], 12)
            self.assertEqual(payload["top_degree_centrality"][0]["node"], "MAPK1")


if __name__ == "__main__":
    unittest.main()
