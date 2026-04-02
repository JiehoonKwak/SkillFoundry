from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHEM_PYTHON = ROOT / "slurm" / "envs" / "chem-tools" / "bin" / "python"


class Phase40FrontierClosureSkillTests(unittest.TestCase):
    def test_openalex_citation_chain_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openalex_citation_chain.json"
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/openalex-citation-chain-starter/scripts/run_openalex_citation_chain.py",
                    "--work-id",
                    "10.1038/nature12373",
                    "--limit",
                    "3",
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
            self.assertIn(payload["result_origin"], {"live_api", "asset_fallback"})
            self.assertEqual(payload["requested_limit"], 3)
            self.assertGreaterEqual(len(payload["top_citing_works"]), 1)

    def test_reactome_pathway_analysis_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "reactome_pathway_analysis.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/reactome-pathway-analysis-starter/scripts/run_reactome_pathway_analysis.py",
                    "--identifiers",
                    "BRCA1,TP53,EGFR",
                    "--top-n",
                    "5",
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
            self.assertGreaterEqual(payload["significant_pathway_count"], 1)
            self.assertGreaterEqual(len(payload["top_pathways"]), 1)

    def test_networkx_propagation_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "networkx_propagation.json"
            subprocess.run(
                [
                    "python3",
                    "skills/systems-biology/networkx-network-propagation-starter/scripts/run_networkx_network_propagation.py",
                    "--input",
                    "skills/systems-biology/networkx-network-propagation-starter/examples/toy_network.tsv",
                    "--seeds",
                    "skills/systems-biology/networkx-network-propagation-starter/examples/toy_seeds.txt",
                    "--top-k",
                    "5",
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
            self.assertEqual(payload["top_non_seed_nodes"][0]["node"], "GRB2")

    def test_rdkit_conformer_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "rdkit_conformers.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/rdkit-conformer-generation-starter/scripts/run_rdkit_conformer_generation.py",
                    "--input",
                    "skills/computational-chemistry-and-molecular-simulation/rdkit-conformer-generation-starter/examples/molecules.tsv",
                    "--num-confs",
                    "4",
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
            self.assertEqual(payload["molecules"][0]["conformer_count"], 4)
            self.assertEqual(payload["molecules"][1]["conformer_count"], 4)

    def test_openmm_forcefield_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "openmm_forcefield.json"
            subprocess.run(
                [
                    str(CHEM_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/openmm-forcefield-assignment-starter/scripts/run_openmm_forcefield_assignment.py",
                    "--input",
                    "skills/computational-chemistry-and-molecular-simulation/openmm-forcefield-assignment-starter/examples/two_waters.pdb",
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
            self.assertEqual(payload["particle_count"], 6)
            self.assertIn("NonbondedForce", payload["force_classes"])


if __name__ == "__main__":
    unittest.main()
