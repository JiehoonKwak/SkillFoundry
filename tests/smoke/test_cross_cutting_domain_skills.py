from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STATISTICS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class CrossCuttingDomainSkillSmokeTests(unittest.TestCase):
    def test_skill_registry_router_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "router.json"
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-agents-and-automation/skill-registry-router-starter/scripts/route_skill_query.py",
                    "--query",
                    "single-cell marker ranking",
                    "--top-k",
                    "3",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["matches"][0]["slug"], "scanpy-ranked-genes-starter")

    def test_matplotlib_publication_plot_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            png_out = Path(tmp_dir) / "plot.png"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(STATISTICS_PYTHON),
                    "skills/visualization-and-reporting/matplotlib-publication-plot-starter/scripts/render_publication_plot.py",
                    "--input",
                    "skills/visualization-and-reporting/matplotlib-publication-plot-starter/examples/toy_measurements.tsv",
                    "--png-out",
                    str(png_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertTrue(png_out.exists())
            self.assertGreater(payload["r_squared"], 0.98)

    def test_registry_link_audit_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "audit.json"
            subprocess.run(
                [
                    "python3",
                    "skills/meta-maintenance/registry-link-audit-starter/scripts/audit_registry_links.py",
                    "--resource-id",
                    "matplotlib-docs",
                    "--resource-id",
                    "lychee-docs",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["failing_count"], 0)


if __name__ == "__main__":
    unittest.main()
