from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "systems-biology" / "clusterprofiler-custom-enrichment" / "scripts" / "run_clusterprofiler_custom_enrichment.R"
R_LIB = os.environ.get("BIOC_SKILL_R_LIB", "")


class ClusterProfilerCustomEnrichmentTests(unittest.TestCase):
    def test_describe_toy_reports_expected_counts(self) -> None:
        completed = subprocess.run(
            ["Rscript", str(SCRIPT), "--describe-toy"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["query_gene_count"], 4)
        self.assertEqual(payload["term_count"], 3)
        self.assertEqual(payload["term_names"]["TERM_UP"], "Positive signal term")

    def test_clusterprofiler_run_with_optional_local_library(self) -> None:
        if not R_LIB:
            self.skipTest("Set BIOC_SKILL_R_LIB to a library containing clusterProfiler to run the package-backed integration test.")
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "clusterprofiler.json"
            subprocess.run(
                ["Rscript", str(SCRIPT), "--lib-dir", R_LIB, "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["package"], "clusterProfiler")
            self.assertGreaterEqual(payload["result_count"], 1)
            self.assertEqual(payload["results"][0]["term"], "TERM_UP")
            self.assertIn("GENE001", payload["results"][0]["gene_ids"])


if __name__ == "__main__":
    unittest.main()
