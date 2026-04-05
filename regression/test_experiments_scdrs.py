from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_ROOT = ROOT / "experiments" / "scDRS"
SKILL_ROOT = EXPERIMENT_ROOT / "scdrs_gwas_to_singlecell_skill"


class ScDRSExperimentContractTests(unittest.TestCase):
    def test_required_experiment_files_exist(self) -> None:
        required_paths = [
            EXPERIMENT_ROOT / "README.md",
            SKILL_ROOT / "SKILL.md",
            SKILL_ROOT / "refs.md",
            SKILL_ROOT / "resource_inventory.md",
            SKILL_ROOT / "deliverables_checklist.md",
            SKILL_ROOT / "templates" / "scdrs_execution_checklist.md",
        ]
        for path in required_paths:
            self.assertTrue(path.is_file(), f"Missing required experiment file: {path}")

    def test_skill_mentions_required_workflow_stages_and_failure_modes(self) -> None:
        skill_text = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
        required_phrases = [
            "scdrs munge-gs",
            "scdrs compute-score",
            "scdrs perform-downstream",
            "positive-control",
            "If `magma` is missing",
            "If `scdrs` is missing",
            "If the GWAS format is incompatible",
            "If `adata.obs` lacks needed annotations",
            "10-kb window around the gene body",
            "top 1,000 genes",
            "HLA/MHC region",
            "--flag-filter-data",
            "dash and underscore spellings as equivalent",
            "Benjamini-Hochberg correction",
            "cell-type associations at `FDR < 0.05`",
            "individual significant cells at `FDR < 0.1`",
            "continuous numeric annotations in `adata.obs`",
            "fewer than 100 cells as exploratory or low-power",
        ]
        for phrase in required_phrases:
            self.assertIn(phrase, skill_text)

    def test_reference_inventory_covers_requested_canonical_sources(self) -> None:
        refs_text = (SKILL_ROOT / "refs.md").read_text(encoding="utf-8")
        inventory_text = (SKILL_ROOT / "resource_inventory.md").read_text(encoding="utf-8")
        expected_refs = [
            "martinjzhang.github.io/scDRS",
            "github.com/martinjzhang/scDRS",
            "github.com/martinjzhang/scDRS/issues/2",
            "pmc.ncbi.nlm.nih.gov/articles/PMC9891382/",
            "nature.com/articles/s41588-022-01167-z",
            "ctg.cncr.nl/software/magma",
            "ebi.ac.uk/gwas",
            "github.com/sulab-wmu/scPagwas",
        ]
        for expected in expected_refs:
            self.assertIn(expected, refs_text)
            self.assertIn(expected, inventory_text)

    def test_checklists_define_required_output_contract(self) -> None:
        deliverables_text = (SKILL_ROOT / "deliverables_checklist.md").read_text(encoding="utf-8")
        execution_text = (SKILL_ROOT / "templates" / "scdrs_execution_checklist.md").read_text(
            encoding="utf-8"
        )
        expected_items = [
            "run_config.yaml",
            "trait.genes.out",
            "trait.gs",
            "adata.prepared.h5ad",
            "covariates.cov",
            "trait.score.gz",
            "trait.full_score.gz",
            "positive_control_report.md",
            "final_report.md",
        ]
        for item in expected_items:
            self.assertIn(item, deliverables_text)
            self.assertIn(item, execution_text)

        expected_contract_phrases = [
            "HLA/MHC region",
            "10 kb",
            "top `1,000` genes",
            "--flag-filter-data",
            "continuous numeric",
            "<100` cells",
            "FDR < 0.05",
            "FDR < 0.1",
        ]
        for phrase in expected_contract_phrases:
            self.assertIn(phrase, deliverables_text)
            self.assertIn(phrase, execution_text)


if __name__ == "__main__":
    unittest.main()
