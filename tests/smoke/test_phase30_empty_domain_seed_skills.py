from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCIPY_PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
NEURO_PYTHON = ROOT / "slurm" / "envs" / "neuro" / "bin" / "python"
ASTRO_PYTHON = ROOT / "slurm" / "envs" / "astronomy" / "bin" / "python"
PLANT_PYTHON = ROOT / "slurm" / "envs" / "plant-science" / "bin" / "python"
AUTOMATION_PYTHON = ROOT / "slurm" / "envs" / "automation" / "bin" / "python"


class Phase30EmptyDomainSeedSkillSmokeTests(unittest.TestCase):
    def test_scipy_ode_simulation_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "ode.json"
            subprocess.run(
                [
                    str(SCIPY_PYTHON),
                    "skills/scientific-computing-and-numerical-methods/scipy-ode-simulation-starter/scripts/run_scipy_ode_simulation.py",
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
            self.assertEqual(payload["model"], "lotka_volterra")
            self.assertAlmostEqual(payload["final_state"]["predator"], 2.270535, places=6)

    def test_gbif_species_occurrence_search_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "gbif.json"
            subprocess.run(
                [
                    "python3",
                    "skills/ecology-evolution-and-biodiversity/gbif-species-occurrence-search-starter/scripts/run_gbif_species_occurrence_search.py",
                    "--scientific-name",
                    "Puma concolor",
                    "--country",
                    "US",
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
            self.assertTrue(payload["matched_scientific_name"].startswith("Puma concolor"))
            self.assertGreaterEqual(payload["occurrence_count"], 1)

    def test_nibabel_nifti_summary_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "nifti.json"
            nifti_path = Path(tmp_dir) / "toy.nii.gz"
            subprocess.run(
                [
                    str(NEURO_PYTHON),
                    "skills/neuroscience-and-neuroimaging/nibabel-nifti-summary-starter/scripts/run_nibabel_nifti_summary.py",
                    "--nifti-out",
                    str(nifti_path),
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
            self.assertEqual(payload["shape"], [4, 4, 3, 2])
            self.assertEqual(payload["time_units"], "sec")

    def test_astropy_fits_image_summary_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fits.json"
            fits_path = Path(tmp_dir) / "toy.fits"
            subprocess.run(
                [
                    str(ASTRO_PYTHON),
                    "skills/physics-and-astronomy/astropy-fits-image-summary-starter/scripts/run_astropy_fits_image_summary.py",
                    "--fits-out",
                    str(fits_path),
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
            self.assertEqual(payload["filter"], "r")
            self.assertEqual(payload["shape"], [5, 5])

    def test_plantcv_plant_phenotyping_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "plantcv.json"
            subprocess.run(
                [
                    str(PLANT_PYTHON),
                    "skills/agriculture-food-and-plant-science/plantcv-plant-phenotyping-starter/scripts/run_plantcv_plant_phenotyping.py",
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
            self.assertEqual(payload["foreground_pixel_count"], 768)
            self.assertEqual(payload["bbox_width"], 24)

    def test_opentrons_liquid_handling_protocol_starter(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "opentrons.json"
            protocol_path = Path(tmp_dir) / "protocol.py"
            subprocess.run(
                [
                    str(AUTOMATION_PYTHON),
                    "skills/robotics-lab-automation-and-scientific-instrumentation/opentrons-liquid-handling-protocol-starter/scripts/run_opentrons_liquid_handling_protocol.py",
                    "--protocol-out",
                    str(protocol_path),
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
            self.assertEqual(payload["command_count"], 5)
            self.assertIn("Dropping tip", payload["last_command"])


if __name__ == "__main__":
    unittest.main()
