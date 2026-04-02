# References

This skill now covers a concrete local `gwas` summary-statistics QC path and points to the canonical resources needed once the file is clean enough for downstream analysis.

- `plink2-docs`: [PLINK 2.0 Docs](https://www.cog-genomics.org/plink/2.0/)
  Use for association reruns, clumping, variant filtering, and follow-up locus work.
- `gwas-meta-analysis-qc-protocol`: [Quality control and conduct of genome-wide association meta-analyses](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083217/)
  Canonical consortium-oriented QC protocol covering study-level, meta-level, and output-level checks.
- `ldsc-repo`: [LDSC](https://github.com/bulik/ldsc)
  Use after formatting and harmonization for LD score regression, heritability, and genetic correlation.
- `fuma-docs`: [FUMA GWAS](https://fuma.ctglab.nl/)
  Use after defining lead loci for SNP2GENE and GENE2FUNC interpretation.
- `gwaslab-docs`: [GWASLab](https://cloufield.github.io/gwaslab/)
  Use when the starter flags format drift, palindromic variants, allele harmonization needs, or export to LDSC/FUMA/GWAS-SSF.
- `gwas-catalog-summary-statistics-docs`: [GWAS Catalog Summary Statistics](https://www.ebi.ac.uk/gwas/summary-statistics/docs/)
  Use when you need a canonical submission or interchange format for cleaned output.

Suggested progression:

1. Run the local starter on raw or consortium-delivered summary statistics.
2. Resolve hard failures and harmonization warnings with GWASLab or upstream study owners.
3. Apply clumping or conditional follow-up with PLINK 2.0.
4. Estimate heritability or genetic correlation with LDSC when ancestry-matched LD references are available.
5. Prioritize loci and genes with FUMA or another interpretation layer once lead signals are stable.
