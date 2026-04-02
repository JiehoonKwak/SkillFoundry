from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = ROOT / "skills" / "scientific-knowledge" / "openalex-citation-chain-starter" / "scripts" / "run_openalex_citation_chain.py"
SPEC = importlib.util.spec_from_file_location("openalex_citation_chain", SCRIPT_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class OpenAlexCitationChainTests(unittest.TestCase):
    def test_normalize_work_id_rejects_empty_input(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.normalize_work_id("  ")

    def test_asset_fallback_is_available_for_canonical_target(self) -> None:
        payload = MODULE.load_asset_fallback("10.1038/nature12373", "forced test fallback")
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["result_origin"], "asset_fallback")
        self.assertEqual(payload["fallback_reason"], "forced test fallback")

    def test_live_query_writes_citation_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "citation_chain.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT_PATH),
                    "--work-id",
                    "10.1038/nature12373",
                    "--limit",
                    "2",
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
            self.assertEqual(payload["requested_limit"], 2)
            self.assertGreaterEqual(payload["seed_work"]["cited_by_count"], 1)
            self.assertGreaterEqual(len(payload["top_referenced_works"]), 1)
            self.assertGreaterEqual(len(payload["top_citing_works"]), 1)


if __name__ == "__main__":
    unittest.main()
