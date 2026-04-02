from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "umap-dimensionality-reduction-starter" / "scripts" / "run_umap_dimensionality_reduction.py"
EXAMPLE = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "umap-dimensionality-reduction-starter" / "examples" / "toy_embedding_input.tsv"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"

SPEC = importlib.util.spec_from_file_location("umap_embedding", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class UmapDimensionalityReductionTests(unittest.TestCase):
    def test_load_embedding_input(self) -> None:
        sample_ids, labels, matrix = MODULE.load_embedding_input(EXAMPLE)
        self.assertEqual(len(sample_ids), 6)
        self.assertEqual(matrix.shape, (6, 4))
        self.assertEqual(labels.count("alpha"), 3)

    def test_run_embedding(self) -> None:
        payload = MODULE.run_embedding(EXAMPLE)
        self.assertEqual(payload["sample_count"], 6)
        self.assertEqual(payload["labels"], ["alpha", "beta"])
        self.assertGreater(payload["centroid_distance"], 1.0)

    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(STATS_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(len(payload["points"]), 6)


if __name__ == "__main__":
    unittest.main()
