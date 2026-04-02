from __future__ import annotations

import importlib.util
import io
import json
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = ROOT / "skills/genomics/ensembl-gene-lookup/scripts/lookup_gene.py"


def load_module():
    spec = importlib.util.spec_from_file_location("ensembl_lookup_gene", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class EnsemblLookupFallbackTests(unittest.TestCase):
    def test_cached_asset_is_used_when_live_lookup_fails(self) -> None:
        module = load_module()
        stdout = io.StringIO()
        argv = [
            "lookup_gene.py",
            "--symbol",
            "BRCA2",
            "--species",
            "homo_sapiens",
        ]
        with mock.patch.object(module, "lookup_gene_by_symbol", side_effect=RuntimeError("timed out")):
            with mock.patch.object(sys, "argv", argv):
                with redirect_stdout(stdout):
                    exit_code = module.main()

        payload = json.loads(stdout.getvalue())
        self.assertEqual(exit_code, 0)
        self.assertEqual(payload["lookup"]["display_name"], "BRCA2")
        self.assertEqual(payload["source_mode"], "asset_fallback")
        self.assertEqual(payload["xrefs_source"], "asset")
        self.assertIn("timed out", payload["fallback_reason"])


if __name__ == "__main__":
    unittest.main()
