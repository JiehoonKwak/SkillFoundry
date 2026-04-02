#!/usr/bin/env python3
"""Tests for the PRIDE project search helper."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = SKILL_ROOT / "scripts" / "search_pride_projects.py"


def load_script_module():
    spec = importlib.util.spec_from_file_location("pride_project_search", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class PrideProjectSearchTests(unittest.TestCase):
    def test_normalize_keyword_requires_non_empty_value(self) -> None:
        self.assertEqual(MODULE.normalize_keyword(" phosphoproteomics "), "phosphoproteomics")
        with self.assertRaises(ValueError):
            MODULE.normalize_keyword("   ")

    def test_build_summary_extracts_compact_fields(self) -> None:
        payload = [
            {
                "accession": "PXD000001",
                "title": "Example proteomics project",
                "projectDescription": "Compact project description",
                "submissionType": "COMPLETE",
                "publicationDate": "2026-01-01",
                "submissionDate": "2025-12-20",
                "updatedDate": "2026-01-02",
                "organisms": ["Homo sapiens (human)"],
                "keywords": ["Proteomics", "Mass spectrometry"],
                "experimentTypes": ["Bottom-up proteomics"],
                "references": ["PMID:12345678"],
            }
        ]

        summary = MODULE.build_summary(payload, "proteomics", 1, 0)

        self.assertEqual(summary["project_count"], 1)
        self.assertEqual(summary["projects"][0]["accession"], "PXD000001")
        self.assertEqual(summary["projects"][0]["submission_type"], "COMPLETE")
        self.assertEqual(summary["projects"][0]["organisms"], ["Homo sapiens (human)"])
        self.assertEqual(summary["projects"][0]["references"], ["PMID:12345678"])
        self.assertIn("keyword=proteomics", summary["query_url"])

    def test_live_search_writes_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "pride-projects.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--keyword",
                    "phosphoproteomics",
                    "--page-size",
                    "1",
                    "--out",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )

            payload = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["keyword"], "phosphoproteomics")
            self.assertEqual(payload["project_count"], 1)
            self.assertTrue(payload["projects"][0]["accession"].startswith("PXD"))
            self.assertTrue(payload["projects"][0]["title"])
            self.assertIn("publication_date", payload["projects"][0])
            self.assertIn("submission_date", payload["projects"][0])
            self.assertIn("references", payload["projects"][0])
            self.assertTrue(payload["projects"][0]["project_url"].endswith(payload["projects"][0]["accession"]))


if __name__ == "__main__":
    unittest.main()
