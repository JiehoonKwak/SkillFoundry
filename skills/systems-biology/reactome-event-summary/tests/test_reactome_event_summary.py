#!/usr/bin/env python3
"""Tests for the Reactome event summary helper."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = SKILL_ROOT / "scripts" / "fetch_reactome_event_summary.py"


def load_script_module():
    spec = importlib.util.spec_from_file_location("reactome_event_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class ReactomeEventSummaryTests(unittest.TestCase):
    def test_build_summary_extracts_expected_fields(self) -> None:
        payload = {
            "stId": "R-HSA-TEST",
            "stIdVersion": "R-HSA-TEST.1",
            "dbId": 123,
            "displayName": "Example pathway",
            "schemaClass": "Pathway",
            "speciesName": "Homo sapiens",
            "reviewStatus": "reviewed",
            "isInferred": False,
            "isInDisease": False,
            "compartment": [{"displayName": "cytosol"}],
            "input": [{"stId": "R-HSA-IN", "displayName": "Input molecule", "schemaClass": "EntityWithAccessionedSequence"}],
            "output": [{"stId": "R-HSA-OUT", "displayName": "Output molecule", "schemaClass": "SimpleEntity"}],
            "literatureReference": [{"pubMedIdentifier": 12345678, "displayName": "Example paper"}],
            "summation": [{"text": "Example summary text."}],
        }

        summary = MODULE.build_summary(payload)

        self.assertEqual(summary["stable_id"], "R-HSA-TEST")
        self.assertEqual(summary["display_name"], "Example pathway")
        self.assertEqual(summary["schema_class"], "Pathway")
        self.assertEqual(summary["species"], "Homo sapiens")
        self.assertEqual(summary["review_status"], "reviewed")
        self.assertEqual(summary["compartments"], ["cytosol"])
        self.assertEqual(summary["inputs"][0]["stable_id"], "R-HSA-IN")
        self.assertEqual(summary["outputs"][0]["display_name"], "Output molecule")
        self.assertEqual(summary["summary_text"], "Example summary text.")

    def test_live_fetch_writes_json_for_known_stable_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--stable-id",
                    "R-HSA-141409",
                    "--out",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            payload = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["stable_id"], "R-HSA-141409")
            self.assertEqual(payload["display_name"], "Mad1 binds kinetochore")
            self.assertIn(payload["schema_class"], {"Reaction", "Pathway"})
            self.assertEqual(payload["species"], "Homo sapiens")
            self.assertTrue(payload["review_status"])
            self.assertTrue(payload["source_url"].endswith("/R-HSA-141409"))

    def test_asset_fallback_is_available_for_canonical_stable_id(self) -> None:
        payload = MODULE.load_asset_fallback("R-HSA-141409", "forced test fallback")
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["result_origin"], "asset_fallback")
        self.assertEqual(payload["fallback_reason"], "forced test fallback")
        self.assertEqual(payload["stable_id"], "R-HSA-141409")
        self.assertEqual(payload["display_name"], "Mad1 binds kinetochore")


if __name__ == "__main__":
    unittest.main()
