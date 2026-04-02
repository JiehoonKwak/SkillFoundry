from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = ROOT / "skills" / "systems-biology" / "reactome-pathway-analysis-starter" / "scripts" / "run_reactome_pathway_analysis.py"
SPEC = importlib.util.spec_from_file_location("reactome_pathway_analysis", SCRIPT_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ReactomePathwayAnalysisTests(unittest.TestCase):
    def test_normalize_identifiers_deduplicates(self) -> None:
        identifiers = MODULE.normalize_identifiers(raw=" BRCA1,TP53,BRCA1,EGFR ")
        self.assertEqual(identifiers, ["BRCA1", "TP53", "EGFR"])

    def test_asset_fallback_is_available_for_canonical_input(self) -> None:
        payload = MODULE.load_asset_fallback(["BRCA1", "TP53", "EGFR"], "forced test fallback")
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["result_origin"], "asset_fallback")

    def test_live_query_writes_ranked_pathway_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "reactome_pathway_analysis.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT_PATH),
                    "--identifiers",
                    "BRCA1,TP53,EGFR",
                    "--top-n",
                    "5",
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
            self.assertIn(payload["result_origin"], {"live_api", "asset_fallback"})
            self.assertEqual(payload["identifiers"], ["BRCA1", "TP53", "EGFR"])
            self.assertGreaterEqual(len(payload["top_pathways"]), 1)
            self.assertGreaterEqual(payload["significant_pathway_count"], 1)
            self.assertIsNotNone(payload["best_fdr"])


if __name__ == "__main__":
    unittest.main()
