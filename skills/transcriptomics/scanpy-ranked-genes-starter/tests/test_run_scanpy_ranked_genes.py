from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCANPY_PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "scanpy-ranked-genes-starter"
    / "scripts"
    / "run_scanpy_ranked_genes.py"
)
COUNT_MATRIX = ROOT / "skills" / "transcriptomics" / "scanpy-ranked-genes-starter" / "examples" / "toy_counts.tsv"
GROUPS = ROOT / "skills" / "transcriptomics" / "scanpy-ranked-genes-starter" / "examples" / "toy_groups.tsv"


def load_script_module():
    spec = importlib.util.spec_from_file_location("scanpy_ranked_genes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class ScanpyRankedGenesTests(unittest.TestCase):
    def test_load_groups_returns_expected_mapping(self) -> None:
        groups = MODULE.load_groups(GROUPS)
        self.assertEqual(groups["cell_a1"], "group_a")
        self.assertEqual(groups["cell_b3"], "group_b")

    def test_runtime_ranked_genes_summary(self) -> None:
        self.assertTrue(SCANPY_PYTHON.exists(), f"Missing Scanpy interpreter: {SCANPY_PYTHON}")
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "ranked-genes.json"
            subprocess.run(
                [
                    str(SCANPY_PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(COUNT_MATRIX),
                    "--groups",
                    str(GROUPS),
                    "--top-n",
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
            self.assertEqual(payload["cells"], 6)
            self.assertEqual(payload["genes"], 6)
            self.assertEqual(payload["groups"], ["group_a", "group_b"])
            self.assertEqual(payload["top_markers_by_group"]["group_a"][0]["gene"], "GeneA")
            self.assertEqual(payload["top_markers_by_group"]["group_b"][0]["gene"], "GeneC")


if __name__ == "__main__":
    unittest.main()
