from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class WorkflowRuntimeSkillSmokeTests(unittest.TestCase):
    def test_nextflow_hello_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            summary_path = tmp_path / "nextflow-summary.json"
            out_dir = tmp_path / "outputs"
            work_dir = tmp_path / "work"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/nextflow-hello-workflow/scripts/run_nextflow_hello.py",
                    "--out-dir",
                    str(out_dir),
                    "--work-dir",
                    str(work_dir),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["executor"], "local")
            self.assertEqual(len(payload["files"]), 4)
            self.assertTrue(any("Version:" in line for line in payload["nextflow_info"]))
            self.assertEqual(payload["slurm_jobs"], [])
            self.assertEqual(sorted(item["name"] for item in payload["files"]), ["Bonjour.txt", "Ciao.txt", "Hello.txt", "Hola.txt"])

    def test_nf_core_pipeline_list(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "nfcore-summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/nf-core-pipeline-list/scripts/list_nfcore_pipelines.py",
                    "--sort",
                    "pulled",
                    "--limit",
                    "3",
                    "--out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertGreaterEqual(payload["counts"]["remote_workflows"], 1)
            self.assertEqual(len(payload["remote_workflows"]), 3)
            self.assertIn("local_workflows", payload["counts"])


if __name__ == "__main__":
    unittest.main()
