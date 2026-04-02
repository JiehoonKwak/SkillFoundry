from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "visualization-and-reporting" / "mkdocs-summary-catalog-starter" / "scripts" / "build_mkdocs_summary_catalog.py"
EXAMPLE = ROOT / "skills" / "visualization-and-reporting" / "mkdocs-summary-catalog-starter" / "examples" / "toy_catalog.json"
REPORT_PYTHON = ROOT / "slurm" / "envs" / "reporting" / "bin" / "python"

SPEC = importlib.util.spec_from_file_location("mkdocs_catalog", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class MkdocsSummaryCatalogTests(unittest.TestCase):
    def test_run_build(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = MODULE.run_build(EXAMPLE, Path(tmp_dir) / "workspace")
            self.assertGreaterEqual(payload["html_file_count"], 3)
            self.assertIn("index.html", payload["html_files"])

    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workspace = Path(tmp_dir) / "workspace"
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(REPORT_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--workspace", str(workspace), "--summary-out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertIn("index.html", payload["html_files"])


if __name__ == "__main__":
    unittest.main()
