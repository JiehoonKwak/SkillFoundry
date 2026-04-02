from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA_PYTHON = ROOT / "slurm" / "envs" / "data-tools" / "bin" / "python"
GENOMICS_PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"
INSTRUMENTATION_PYTHON = ROOT / "slurm" / "envs" / "instrumentation" / "bin" / "python"


class Phase37ResourceConversionSmokeTests(unittest.TestCase):
    def test_frictionless_tabular_validation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    "skills/data-acquisition-and-dataset-handling/frictionless-tabular-validation-starter/scripts/run_frictionless_tabular_validation.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/frictionless-tabular-validation-starter/examples/toy_people_valid.csv",
                    "--schema",
                    "skills/data-acquisition-and-dataset-handling/frictionless-tabular-validation-starter/examples/toy_people_schema.json",
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
            self.assertTrue(payload["valid"])
            self.assertEqual(payload["row_count"], 3)

    def test_rocrate_metadata_bundle_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            out_path = tmp_root / "summary.json"
            crate_dir = tmp_root / "toy_bundle"
            subprocess.run(
                [
                    str(DATA_PYTHON),
                    "skills/data-acquisition-and-dataset-handling/rocrate-metadata-bundle-starter/scripts/build_rocrate_metadata_bundle.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/rocrate-metadata-bundle-starter/examples/toy_measurements.csv",
                    "--crate-dir",
                    str(crate_dir),
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertTrue(payload["has_metadata_file"])
            self.assertEqual(payload["data_entity_count"], 1)
            self.assertEqual(payload["bundled_files"], ["toy_measurements.csv"])

    def test_minimap2_read_mapping_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            out_path = tmp_root / "summary.json"
            bam_out = tmp_root / "toy_reads.bam"
            subprocess.run(
                [
                    str(GENOMICS_PYTHON),
                    "skills/genomics/minimap2-read-mapping-starter/scripts/run_minimap2_read_mapping.py",
                    "--reference",
                    "skills/genomics/minimap2-read-mapping-starter/examples/toy_reference.fa",
                    "--reads",
                    "skills/genomics/minimap2-read-mapping-starter/examples/toy_reads.fastq",
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
            self.assertEqual(payload["mapped_count"], 2)
            self.assertEqual(payload["unmapped_count"], 1)
            self.assertEqual(payload["reference_names"], ["chrToy"])

    def test_qcodes_parameter_sweep_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            out_path = tmp_root / "summary.json"
            db_out = tmp_root / "toy_sweep.db"
            subprocess.run(
                [
                    str(INSTRUMENTATION_PYTHON),
                    "skills/robotics-lab-automation-and-scientific-instrumentation/qcodes-parameter-sweep-starter/scripts/run_qcodes_parameter_sweep.py",
                    "--setpoints",
                    "skills/robotics-lab-automation-and-scientific-instrumentation/qcodes-parameter-sweep-starter/examples/toy_setpoints.tsv",
                    "--db-out",
                    str(db_out),
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
            self.assertEqual(payload["point_count"], 4)
            self.assertAlmostEqual(payload["response_mean"], 1.15, places=6)
            self.assertEqual(payload["run_id"], 1)


if __name__ == "__main__":
    unittest.main()
