from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "systems-biology" / "fgsea-preranked-enrichment" / "scripts" / "run_fgsea_preranked.R"
RSCRIPT = ROOT / "slurm" / "envs" / "bioconductor" / "bin" / "Rscript"
R_LIB = os.environ.get("BIOC_SKILL_R_LIB", "")


class FgseaPrerankedEnrichmentTests(unittest.TestCase):
    def test_describe_toy_reports_expected_sizes(self) -> None:
        completed = subprocess.run(
            [str(RSCRIPT), str(SCRIPT), "--describe-toy"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["gene_count"], 12)
        self.assertEqual(payload["pathway_count"], 3)
        self.assertEqual(payload["strongest_positive_gene"], "GENE001")
        self.assertEqual(payload["pathway_sizes"]["positive_signal"], 4)

    def test_fgsea_run_with_optional_local_library(self) -> None:
        if not R_LIB:
            self.skipTest("Set BIOC_SKILL_R_LIB to a library containing fgsea to run the package-backed integration test.")
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "fgsea.json"
            completed = subprocess.run(
                [str(RSCRIPT), str(SCRIPT), "--lib-dir", R_LIB, "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            self.assertEqual(completed.stderr, "")
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["package"], "fgsea")
            self.assertGreaterEqual(payload["result_count"], 1)
            pathways = [row["pathway"] for row in payload["results"]]
            self.assertIn("positive_signal", pathways)


if __name__ == "__main__":
    unittest.main()
