from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENOMICS_PYTHON = ROOT / "slurm" / "envs" / "genomics" / "bin" / "python"
NEURO_PYTHON = ROOT / "slurm" / "envs" / "neuro" / "bin" / "python"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class Phase33FrontierSkillConversionSmokeTests(unittest.TestCase):
    def test_fastqc_multiqc_read_qc_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fastqc_multiqc_summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/genomics/fastqc-multiqc-read-qc-starter/scripts/run_fastqc_multiqc_read_qc.py",
                    "--input",
                    "skills/genomics/fastqc-multiqc-read-qc-starter/examples/toy_reads.fastq",
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["total_sequences"], 4)
            self.assertEqual(payload["gc_percent"], 50)
            self.assertEqual(payload["multiqc_sample_count"], 1)

    def test_macs3_peak_calling_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "macs3_peak_summary.json"
            subprocess.run(
                [
                    str(GENOMICS_PYTHON),
                    "skills/epigenomics-and-chromatin/macs3-peak-calling-starter/scripts/run_macs3_peak_calling.py",
                    "--treatment",
                    "skills/epigenomics-and-chromatin/macs3-peak-calling-starter/examples/toy_treatment.bed",
                    "--summary-out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["peak_count"], 1)
            self.assertEqual(payload["top_peak"]["name"], "toy_peak_1")

    def test_fairlearn_bias_audit_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fairlearn_summary.json"
            subprocess.run(
                [
                    str(STATS_PYTHON),
                    "skills/clinical-biomedical-data-science/fairlearn-bias-audit-starter/scripts/run_fairlearn_bias_audit.py",
                    "--input",
                    "skills/clinical-biomedical-data-science/fairlearn-bias-audit-starter/examples/toy_fairness_cohort.tsv",
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
            self.assertEqual(payload["fairness_metrics"]["demographic_parity_difference"], 0.5)
            self.assertTrue(payload["flags"])

    def test_mne_connectivity_graph_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "mne_connectivity_graph.json"
            subprocess.run(
                [
                    str(NEURO_PYTHON),
                    "skills/neuroscience-and-neuroimaging/mne-connectivity-graph-starter/scripts/run_mne_connectivity_graph.py",
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
            self.assertEqual(payload["max_edge"], {"source": "Fz", "target": "Cz", "weight": 0.926824})
            self.assertEqual(len(payload["strong_edges"]), 1)

    def test_uniprot_sequence_feature_annotation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "uniprot_summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/proteomics/uniprot-sequence-feature-annotation-starter/scripts/fetch_uniprot_sequence_feature_summary.py",
                    "--accession",
                    "P04637",
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
            self.assertEqual(payload["accession"], "P04637")
            self.assertGreaterEqual(payload["feature_count"], 1000)
            self.assertEqual(payload["representative_features"][0]["type"], "Chain")


if __name__ == "__main__":
    unittest.main()
