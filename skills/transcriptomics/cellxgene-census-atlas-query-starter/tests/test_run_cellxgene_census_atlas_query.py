from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "census" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "cellxgene-census-atlas-query-starter"
    / "scripts"
    / "run_cellxgene_census_atlas_query.py"
)


class CellxgeneCensusAtlasQueryTests(unittest.TestCase):
    def test_runtime_query_returns_tabula_sapiens_matches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "census_query.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--keyword",
                    "Tabula Sapiens",
                    "--limit",
                    "5",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertGreaterEqual(payload["match_count"], 1)
            self.assertEqual(payload["top_matches"][0]["collection_name"], "Tabula Sapiens")
            self.assertGreater(payload["top_matches"][0]["dataset_total_cell_count"], 1000000)

    def test_rejects_unknown_keyword(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--keyword",
                "definitely-no-such-census-dataset-keyword",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=300,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("No CELLxGENE Census datasets matched keyword", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
