from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "zarr-chunked-array-store-starter"
    / "scripts"
    / "run_zarr_chunked_array_store.py"
)
INPUT = (
    ROOT
    / "skills"
    / "data-acquisition-and-dataset-handling"
    / "zarr-chunked-array-store-starter"
    / "examples"
    / "toy_matrix.tsv"
)


class ZarrChunkedArrayStoreTests(unittest.TestCase):
    def test_writes_chunked_store_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            summary_out = tmp_root / "summary.json"
            store_out = tmp_root / "toy_matrix.zarr"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--store-out",
                    str(store_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["shape"], [6, 4])
            self.assertEqual(payload["chunk_shape"], [2, 2])
            self.assertEqual(payload["matrix_sum"], 300.0)
            self.assertEqual(payload["first_chunk_sum"], 14.0)
            self.assertGreaterEqual(payload["store_file_count"], 2)

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                "missing.tsv",
                "--store-out",
                "scratch/zarr/missing.zarr",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input matrix not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
