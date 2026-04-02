from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run_json_command(*args: str) -> dict:
    completed = subprocess.run(
        ["python3", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return json.loads(completed.stdout)


class ProteomicsAndCheminformaticsApiSkillSmokeTests(unittest.TestCase):
    def test_pride_project_search_returns_projects(self) -> None:
        payload = run_json_command(
            "skills/proteomics/pride-project-search/scripts/search_pride_projects.py",
            "--keyword",
            "phosphoproteomics",
            "--page-size",
            "1",
        )
        self.assertEqual(payload["keyword"], "phosphoproteomics")
        self.assertGreaterEqual(payload["project_count"], 1)
        self.assertTrue(payload["projects"][0]["accession"].startswith("PXD"))

    def test_ebi_proteins_entry_summary_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/proteomics/ebi-proteins-entry-summary/scripts/fetch_protein_summary.py",
            "--accession",
            "P38398",
        )
        self.assertEqual(payload["accession"], "P38398")
        self.assertIn("BRCA1", payload["gene_names"])
        self.assertEqual(payload["organism"], "Homo sapiens")
        self.assertTrue(payload["source_url"].endswith("/P38398"))

    def test_chembl_molecule_search_returns_known_compound(self) -> None:
        payload = run_json_command(
            "skills/drug-discovery-and-cheminformatics/chembl-molecule-search/scripts/search_chembl_molecules.py",
            "--query",
            "imatinib",
            "--limit",
            "1",
        )
        self.assertEqual(payload["query"], "imatinib")
        self.assertEqual(len(payload["molecules"]), 1)
        self.assertEqual(payload["molecules"][0]["chembl_id"], "CHEMBL941")
        self.assertEqual(payload["molecules"][0]["preferred_name"], "IMATINIB")


if __name__ == "__main__":
    unittest.main()
