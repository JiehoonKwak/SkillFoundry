#!/usr/bin/env python3
"""Tests for the CRISPR screen analysis skill."""

import json
import subprocess
import sys
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "scripts" / "run_crispr_screen.py"
EXAMPLE_COUNTS = SKILL_DIR / "examples" / "toy_counts.tsv"


class TestCrisprScreenSkill(unittest.TestCase):
    """Validate the CRISPR screen analysis on toy data."""

    @classmethod
    def setUpClass(cls):
        """Run the analysis once and load results."""
        cls.out_json = Path("/tmp/test_crispr_screen_results.json")
        cls.out_tsv = Path("/tmp/test_crispr_screen_summary.tsv")
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--counts", str(EXAMPLE_COUNTS),
                "--control", "ctrl_rep1,ctrl_rep2",
                "--treatment", "treatment_rep1,treatment_rep2",
                "--json-out", str(cls.out_json),
                "--tsv-out", str(cls.out_tsv),
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Script failed:\n{result.stderr}")
        cls.payload = json.loads(cls.out_json.read_text())
        cls.gene_map = {
            g["gene"]: g for g in cls.payload["gene_summary"]
        }

    def test_exits_successfully(self):
        """Script should produce valid JSON output."""
        self.assertIn("gene_summary", self.payload)

    def test_gene_count(self):
        """Should detect all 10 genes (9 targets + 1 control)."""
        # Includes: TP53, RB1, PTEN, CDKN2A, MYC, KRAS, EGFR, ACTB, GAPDH, NontargetingControl
        self.assertEqual(self.payload["n_genes"], 10)

    def test_guide_count(self):
        """Should detect all 32 guides."""
        self.assertEqual(self.payload["n_guides"], 32)

    def test_tumor_suppressors_depleted(self):
        """TP53, RB1, PTEN, CDKN2A should rank in top 4 for negative selection."""
        tsg_genes = {"TP53", "RB1", "PTEN", "CDKN2A"}
        top4 = {g["gene"] for g in self.payload["gene_summary"][:4]}
        self.assertEqual(top4, tsg_genes)

    def test_oncogenes_enriched(self):
        """MYC, KRAS, EGFR should rank in top 3 for positive selection."""
        for gene_name in ["MYC", "KRAS", "EGFR"]:
            gene = self.gene_map[gene_name]
            self.assertLessEqual(
                gene["pos"]["rank"], 3,
                f"{gene_name} pos rank should be <= 3, got {gene['pos']['rank']}",
            )

    def test_neutral_genes_not_significant(self):
        """ACTB should not be significant in negative selection."""
        actb = self.gene_map["ACTB"]
        self.assertGreater(
            actb["neg"]["p_value"], 0.05,
            f"ACTB neg p-value should be > 0.05, got {actb['neg']['p_value']}",
        )

    def test_log2fc_direction(self):
        """Depleted genes should have negative LFC, enriched should have positive."""
        for gene_name in ["TP53", "RB1", "PTEN", "CDKN2A"]:
            gene = self.gene_map[gene_name]
            self.assertLess(gene["neg"]["lfc"], -1.0, f"{gene_name} LFC should be < -1")

        for gene_name in ["MYC", "KRAS", "EGFR"]:
            gene = self.gene_map[gene_name]
            self.assertGreater(gene["pos"]["lfc"], 1.0, f"{gene_name} LFC should be > 1")

    def test_size_factors_reasonable(self):
        """Size factors should be close to 1.0 for balanced libraries."""
        for sample, sf in self.payload["size_factors"].items():
            self.assertGreater(sf, 0.5, f"Size factor for {sample} too low: {sf}")
            self.assertLess(sf, 2.0, f"Size factor for {sample} too high: {sf}")

    def test_fdr_bounded(self):
        """All FDR values should be between 0 and 1."""
        for gene in self.payload["gene_summary"]:
            self.assertGreaterEqual(gene["neg"]["fdr"], 0.0)
            self.assertLessEqual(gene["neg"]["fdr"], 1.0)
            self.assertGreaterEqual(gene["pos"]["fdr"], 0.0)
            self.assertLessEqual(gene["pos"]["fdr"], 1.0)

    def test_tsv_output_exists(self):
        """TSV gene summary should be written."""
        self.assertTrue(self.out_tsv.exists())
        lines = self.out_tsv.read_text().strip().splitlines()
        # Header + 10 genes
        self.assertEqual(len(lines), 11)

    def test_guide_details_present(self):
        """Each gene should have guide-level details."""
        tp53 = self.gene_map["TP53"]
        self.assertEqual(len(tp53["guides"]), 4)
        for guide in tp53["guides"]:
            self.assertIn("guide", guide)
            self.assertIn("log2fc", guide)
            self.assertIn("p_value", guide)


if __name__ == "__main__":
    unittest.main()
