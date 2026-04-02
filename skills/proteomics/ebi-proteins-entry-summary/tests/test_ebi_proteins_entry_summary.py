#!/usr/bin/env python3
"""Tests for the EBI Proteins entry summary helper."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = SKILL_ROOT / "scripts" / "fetch_protein_entry_summary.py"


def load_module():
    spec = importlib.util.spec_from_file_location("ebi_proteins_entry_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


MODULE = load_module()


class EbiProteinsEntrySummaryTests(unittest.TestCase):
    def test_build_summary_extracts_compact_fields(self) -> None:
        payload = {
            "accession": "PTEST",
            "id": "TEST_HUMAN",
            "protein": {"recommendedName": {"fullName": {"value": "Example protein"}}},
            "gene": [{"name": {"value": "TEST1"}, "synonyms": [{"value": "ALT1"}]}],
            "organism": {
                "taxonomy": 9606,
                "names": [{"type": "scientific", "value": "Homo sapiens"}, {"type": "common", "value": "Human"}],
            },
            "sequence": {"length": 123, "mass": 4567},
            "keywords": [{"value": "DNA-binding"}],
            "comments": [{"type": "FUNCTION", "text": [{"value": "Example function"}]}],
            "features": [{"type": "DOMAIN", "category": "DOMAINS_AND_SITES", "description": "Example domain", "begin": "1", "end": "20"}],
        }

        summary = MODULE.build_summary(payload)

        self.assertEqual(summary["accession"], "PTEST")
        self.assertEqual(summary["recommended_name"], "Example protein")
        self.assertEqual(summary["gene_names"], ["TEST1", "ALT1"])
        self.assertEqual(summary["organism_scientific_name"], "Homo sapiens")
        self.assertEqual(summary["sequence_length"], 123)
        self.assertEqual(summary["keywords"], ["DNA-binding"])
        self.assertEqual(summary["comments"][0]["type"], "FUNCTION")

    def test_live_fetch_writes_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "protein-summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--accession",
                    "P38398",
                    "--out",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )

            payload = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["accession"], "P38398")
            self.assertEqual(payload["organism_scientific_name"], "Homo sapiens")
            self.assertIn("BRCA1", payload["gene_names"])
            self.assertGreaterEqual(payload["sequence_length"], 1800)
            self.assertTrue(payload["recommended_name"])


if __name__ == "__main__":
    unittest.main()
