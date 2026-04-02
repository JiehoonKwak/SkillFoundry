from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"
DEEPCHEM_PYTHON = ROOT / "slurm" / "envs" / "deepchem" / "bin" / "python"
PSI4_PYTHON = ROOT / "slurm" / "envs" / "psi4" / "bin" / "python"
R_SCRIPT = ROOT / "slurm" / "envs" / "bioconductor" / "bin" / "Rscript"


class WorkflowLanguageAndScientificComputeSkillSmokeTests(unittest.TestCase):
    def test_cwl_commandlinetool_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "cwl-summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/cwl-commandlinetool-starter/scripts/run_cwl_hello.py",
                    "--message",
                    "hello from cwl smoke",
                    "--workspace",
                    str(Path(tmp_dir) / "workspace"),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["output_text"], "hello from cwl smoke")

    def test_wdl_task_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "wdl-summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/reproducible-workflows/wdl-task-starter/scripts/run_wdl_hello.py",
                    "--name",
                    "Smoke",
                    "--workspace",
                    str(Path(tmp_dir) / "workspace"),
                    "--summary-out",
                    str(summary_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["greeting_text"], "Hello, Smoke")

    def test_rdkit_molecular_descriptors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "rdkit.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/drug-discovery-and-cheminformatics/rdkit-molecular-descriptors/scripts/compute_rdkit_descriptors.py",
                    "--smiles",
                    "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "--name",
                    "aspirin",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["formula"], "C9H8O4")
            self.assertEqual(payload["hbd"], 1)

    def test_deepchem_circular_featurization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "deepchem.json"
            subprocess.run(
                [
                    str(DEEPCHEM_PYTHON),
                    "skills/drug-discovery-and-cheminformatics/deepchem-circular-featurization/scripts/compute_circular_fingerprints.py",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["featurizer"], "CircularFingerprint")
            self.assertEqual(payload["molecule_count"], 2)
            self.assertEqual(payload["size"], 32)

    def test_openmm_system_minimization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openmm.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/openmm-system-minimization/scripts/run_openmm_minimization.py",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertLess(payload["final_potential_energy_kj_mol"], payload["initial_potential_energy_kj_mol"])

    def test_psi4_single_point_energy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "psi4.json"
            subprocess.run(
                [
                    str(PSI4_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/psi4-single-point-energy/scripts/run_psi4_single_point.py",
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
            self.assertEqual(payload["mode"], "single_point")
            self.assertEqual(payload["method"], "hf")
            self.assertEqual(payload["basis"], "sto-3g")
            self.assertAlmostEqual(payload["energy_hartree"], -74.962991614813, places=9)

    def test_fgsea_preranked_enrichment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fgsea.json"
            subprocess.run(
                [
                    str(R_SCRIPT),
                    "skills/systems-biology/fgsea-preranked-enrichment/scripts/run_fgsea_preranked.R",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["package"], "fgsea")
            self.assertGreaterEqual(payload["result_count"], 1)

    def test_clusterprofiler_custom_enrichment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "clusterprofiler.json"
            subprocess.run(
                [
                    str(R_SCRIPT),
                    "skills/systems-biology/clusterprofiler-custom-enrichment/scripts/run_clusterprofiler_custom_enrichment.R",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["package"], "clusterProfiler")
            self.assertGreaterEqual(payload["result_count"], 1)


if __name__ == "__main__":
    unittest.main()
