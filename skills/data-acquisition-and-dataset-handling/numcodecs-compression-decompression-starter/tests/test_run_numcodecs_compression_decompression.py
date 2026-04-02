from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "numcodecs-compression-decompression-starter" / "scripts" / "run_numcodecs_compression_decompression.py"
EXAMPLE = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "numcodecs-compression-decompression-starter" / "examples" / "toy_matrix.tsv"
DATA_PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"

SPEC = importlib.util.spec_from_file_location("numcodecs_roundtrip", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class NumcodecsCompressionDecompressionTests(unittest.TestCase):
    def test_load_matrix_shape(self) -> None:
        matrix = MODULE.load_matrix(EXAMPLE)
        self.assertEqual(matrix.shape, (3, 4))
        self.assertEqual(str(matrix.dtype), "int32")

    def test_roundtrip_summary(self) -> None:
        payload = MODULE.run_roundtrip(EXAMPLE)
        self.assertTrue(payload["decoded_equal"])
        self.assertEqual(payload["row_sums"], [7, 47, 322])
        self.assertGreater(payload["encoded_nbytes"], 0)

    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLE),
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
            self.assertEqual(payload["codec_name"], "blosc-zstd")
            self.assertTrue(payload["decoded_equal"])


if __name__ == "__main__":
    unittest.main()
