from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    ROOT
    / "skills"
    / "systems-biology"
    / "reactome-identifiers-enrichment"
    / "scripts"
    / "analyze_reactome_identifiers.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("analyze_reactome_identifiers", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class ReactomeIdentifiersEnrichmentTests(unittest.TestCase):
    def test_normalize_identifiers_deduplicates_and_strips(self) -> None:
        module = load_module()
        identifiers = module.normalize_identifiers(raw=" BRCA1,TP53,BRCA1 \n, ")
        self.assertEqual(identifiers, ["BRCA1", "TP53"])

    def test_live_query_writes_expected_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "reactome-enrichment.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/reactome-identifiers-enrichment/scripts/analyze_reactome_identifiers.py",
                    "--identifiers",
                    "BRCA1,TP53",
                    "--page-size",
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
            self.assertEqual(payload["identifiers"], ["BRCA1", "TP53"])
            self.assertGreaterEqual(len(payload["pathways"]), 1)
            self.assertTrue(payload["pathways"][0]["stable_id"].startswith("R-HSA-"))
            self.assertEqual(payload["pathways"][0]["species"], "Homo sapiens")
            self.assertIsNotNone(payload["pathways"][0]["entities_fdr"])


if __name__ == "__main__":
    unittest.main()
