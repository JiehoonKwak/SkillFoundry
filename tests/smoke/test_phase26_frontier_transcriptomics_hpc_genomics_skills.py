from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCANPY_PYTHON = ROOT / "slurm" / "envs" / "scanpy" / "bin" / "python"
CENSUS_PYTHON = ROOT / "slurm" / "envs" / "census" / "bin" / "python"
GENOMICS_PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"


class Phase26TranscriptomicsHpcGenomicsSmokeTests(unittest.TestCase):
    def test_scanpy_combat_batch_correction_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "scanpy_combat.json"
            subprocess.run(
                [
                    str(SCANPY_PYTHON),
                    "skills/transcriptomics/scanpy-combat-batch-correction-starter/scripts/run_scanpy_combat_batch_correction.py",
                    "--counts",
                    "skills/transcriptomics/scanpy-combat-batch-correction-starter/examples/toy_counts.tsv",
                    "--metadata",
                    "skills/transcriptomics/scanpy-combat-batch-correction-starter/examples/toy_metadata.tsv",
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
            self.assertGreater(payload["pre_batch_mean_abs_diff"], payload["post_batch_mean_abs_diff"])
            self.assertLess(payload["batch_gap_ratio"], 0.05)

    def test_cellxgene_census_atlas_query_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "cellxgene_atlas.json"
            subprocess.run(
                [
                    str(CENSUS_PYTHON),
                    "skills/transcriptomics/cellxgene-census-atlas-query-starter/scripts/run_cellxgene_census_atlas_query.py",
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
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertGreaterEqual(payload["match_count"], 1)
            self.assertEqual(payload["top_matches"][0]["collection_name"], "Tabula Sapiens")

    def test_slurm_monitoring_accounting_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "slurm_monitoring.json"
            subprocess.run(
                [
                    "python3",
                    "skills/hpc/slurm-monitoring-accounting-starter/scripts/run_slurm_monitoring_accounting.py",
                    "--partition",
                    "cpu",
                    "--job-name",
                    "phase26-monitoring",
                    "--sleep",
                    "2",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertGreaterEqual(payload["queue_snapshot_count"], 1)
            self.assertEqual(payload["accounting"]["State"], "COMPLETED")
            self.assertEqual(payload["accounting"]["ExitCode"], "0:0")

    def test_pysam_sam_bam_summary_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "pysam_summary.json"
            bam_out = Path(tmp_dir) / "toy_reads.bam"
            subprocess.run(
                [
                    str(GENOMICS_PYTHON),
                    "skills/genomics/pysam-sam-bam-summary-starter/scripts/run_pysam_sam_bam_summary.py",
                    "--input",
                    "skills/genomics/pysam-sam-bam-summary-starter/examples/toy_reads.sam",
                    "--out",
                    str(out_path),
                    "--bam-out",
                    str(bam_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["record_count"], 4)
            self.assertEqual(payload["mapped_count"], 3)
            self.assertEqual(payload["reference_mapped_counts"]["chr1"], 2)
            self.assertTrue(bam_out.exists())
            self.assertTrue((bam_out.parent / (bam_out.name + ".bai")).exists())


if __name__ == "__main__":
    unittest.main()
