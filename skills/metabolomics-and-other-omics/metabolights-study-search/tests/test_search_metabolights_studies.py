from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MODULE_PATH = ROOT / "skills" / "metabolomics-and-other-omics" / "metabolights-study-search" / "scripts" / "search_metabolights_studies.py"
MODULE_SPEC = importlib.util.spec_from_file_location("metabolights_search_module", MODULE_PATH)
assert MODULE_SPEC and MODULE_SPEC.loader
MODULE = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(MODULE)


class MetaboLightsStudySearchTests(unittest.TestCase):
    def test_metabolights_search_returns_at_least_one_study(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "metabolights.json"
            subprocess.run(
                [
                    "python3",
                    "skills/metabolomics-and-other-omics/metabolights-study-search/scripts/search_metabolights_studies.py",
                    "--query",
                    "diabetes",
                    "--rows",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=90,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["query"], "diabetes")
            self.assertEqual(payload["study_count"], 1)
            self.assertGreaterEqual(payload["total_studies"], 1)
            self.assertTrue(payload["studies"][0]["accession"].startswith("MTBLS"))
            self.assertTrue(payload["studies"][0]["title"])

    def test_asset_fallback_is_available_for_canonical_smoke_query(self) -> None:
        payload = MODULE.load_asset_fallback("diabetes", 1, 1, "forced test fallback")
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["result_origin"], "asset_fallback")
        self.assertEqual(payload["fallback_reason"], "forced test fallback")
        self.assertEqual(payload["study_count"], 1)
        self.assertEqual(payload["studies"][0]["accession"], "MTBLS1")


if __name__ == "__main__":
    unittest.main()
