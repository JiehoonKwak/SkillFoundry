from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "pyarrow-format-conversion-starter" / "scripts" / "run_pyarrow_format_conversion.py"
EXAMPLE = ROOT / "skills" / "data-acquisition-and-dataset-handling" / "pyarrow-format-conversion-starter" / "examples" / "toy_matrix.tsv"
DATA_PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"


class PyArrowFormatConversionTests(unittest.TestCase):
    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            parquet_out = Path(tmp_dir) / "toy.parquet"
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLE),
                    "--parquet-out",
                    str(parquet_out),
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
            self.assertEqual(payload["row_count"], 4)
            self.assertEqual(payload["column_count"], 4)
            self.assertTrue(payload["roundtrip_rows_equal"])
            self.assertEqual(payload["schema"][1]["type"], "int64")


if __name__ == "__main__":
    unittest.main()
