from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ECOLOGY_PYTHON = ROOT / "slurm" / "envs" / "ecology" / "bin" / "python"
NEURO_PYTHON = ROOT / "slurm" / "envs" / "neuro" / "bin" / "python"
NUMERICS_PYTHON = ROOT / "slurm" / "envs" / "numerics" / "bin" / "python"


class Phase31FrontierLeafConversionSkillSmokeTests(unittest.TestCase):
    def test_gbif_dataset_search_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "gbif_dataset_search.json"
            subprocess.run(
                [
                    "python3",
                    "skills/ecology-evolution-and-biodiversity/gbif-dataset-search-starter/scripts/run_gbif_dataset_search.py",
                    "--query",
                    "puma",
                    "--limit",
                    "3",
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
            self.assertEqual(payload["query"], "puma")
            self.assertEqual(payload["result_count"], 3)
            self.assertTrue(payload["first_dataset_key"])

    def test_scikitbio_tree_comparison_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "tree_comparison.json"
            subprocess.run(
                [
                    str(ECOLOGY_PYTHON),
                    "skills/ecology-evolution-and-biodiversity/scikitbio-tree-comparison-starter/scripts/run_scikitbio_tree_comparison.py",
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
            self.assertEqual(payload["shared_tip_count"], 4)
            self.assertEqual(payload["robinson_foulds_distance"], 4.0)

    def test_mne_eeg_preprocessing_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "mne_preprocessing.json"
            subprocess.run(
                [
                    str(NEURO_PYTHON),
                    "skills/neuroscience-and-neuroimaging/mne-eeg-preprocessing-starter/scripts/run_mne_eeg_preprocessing.py",
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
            self.assertEqual(payload["channel_names"], ["Fz", "Cz"])
            self.assertEqual(payload["filtered_std_ratio"], [0.931382, 0.977195])

    def test_nilearn_fmri_denoising_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "nilearn_denoising.json"
            subprocess.run(
                [
                    str(NEURO_PYTHON),
                    "skills/neuroscience-and-neuroimaging/nilearn-fmri-denoising-starter/scripts/run_nilearn_fmri_denoising.py",
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
            self.assertEqual(payload["frame_count"], 60)
            self.assertEqual(payload["post_confound_abs_correlation_mean"], 0.0)
            self.assertEqual(payload["cleaned_std_sample"], [1.0, 1.0])

    def test_fipy_diffusion_pde_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fipy_diffusion.json"
            subprocess.run(
                [
                    str(NUMERICS_PYTHON),
                    "skills/scientific-computing-and-numerical-methods/fipy-diffusion-pde-starter/scripts/run_fipy_diffusion_pde.py",
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
            self.assertEqual(payload["cell_count"], 20)
            self.assertEqual(payload["final_mass"], 5.0)
            self.assertEqual(payload["leading_profile"][0], 0.997403)


if __name__ == "__main__":
    unittest.main()
