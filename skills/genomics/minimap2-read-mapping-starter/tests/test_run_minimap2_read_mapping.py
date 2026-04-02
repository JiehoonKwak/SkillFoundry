from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "genomics"
    / "minimap2-read-mapping-starter"
    / "scripts"
    / "run_minimap2_read_mapping.py"
)
EXAMPLES = ROOT / "skills" / "genomics" / "minimap2-read-mapping-starter" / "examples"


class Minimap2ReadMappingTests(unittest.TestCase):
    def test_maps_reads_and_builds_bam(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            bam_out = tmp_root / "toy_reads.bam"
            out_path = tmp_root / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--reference",
                    str(EXAMPLES / "toy_reference.fa"),
                    "--reads",
                    str(EXAMPLES / "toy_reads.fastq"),
                    "--bam-out",
                    str(bam_out),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["read_count"], 3)
            self.assertEqual(payload["mapped_count"], 2)
            self.assertEqual(payload["unmapped_count"], 1)
            self.assertEqual(payload["reference_names"], ["chrToy"])
            self.assertAlmostEqual(payload["mean_mapq"], 26.0, places=6)
            self.assertTrue(bam_out.exists())
            self.assertTrue(Path(str(bam_out) + ".bai").exists())

    def test_missing_reference_fails(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--reference",
                "missing.fa",
                "--reads",
                str(EXAMPLES / "toy_reads.fastq"),
                "--bam-out",
                "scratch/genomics/minimap2/missing.bam",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Reference FASTA not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
