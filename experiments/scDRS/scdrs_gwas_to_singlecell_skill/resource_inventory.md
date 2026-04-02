# Resource Inventory

This inventory focuses on canonical materials beyond the initial user-provided workflow description.

## Core scDRS materials

- scDRS docs home
  URL: https://martinjzhang.github.io/scDRS/
  Why included: official package entry point and install surface.

- scDRS CLI reference
  URL: https://martinjzhang.github.io/scDRS/reference_cli.html
  Why included: portable command surface for the three core stages.

- scDRS file formats
  URL: https://martinjzhang.github.io/scDRS/file_format.html
  Why included: needed to make the skill robust around `.gs`, `.cov`, `.h5ad`, and output contracts.

- scDRS FAQ
  URL: https://martinjzhang.github.io/scDRS/faq.html
  Why included: official statement of the recommended workflow, the top-1,000-gene default, and downstream-analysis expectations.

- scDRS GitHub repository
  URL: https://github.com/martinjzhang/scDRS
  Why included: fallback source for examples, tests, and version-specific behavior.

- scDRS issue #2: custom gene set
  URL: https://github.com/martinjzhang/scDRS/issues/2
  Why included: captures the paper-derived MAGMA setup that this runbook now makes explicit: gene-level p-values, 10-kb SNP-to-gene mapping window, and top 1,000 genes.

## Method papers

- scDRS paper
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC9891382/
  Why included: primary method paper with the clearest accessible statements on MHC exclusion, the MAGMA preprocessing defaults, and reported FDR thresholds.

- scDRS paper publisher version
  URL: https://www.nature.com/articles/s41588-022-01167-z
  Why included: canonical journal landing page for citation metadata and the published article record.

- MAGMA paper
  URL: https://doi.org/10.1371/journal.pcbi.1004219
  Why included: canonical reference for the gene-level GWAS aggregation step feeding scDRS.

- MAGMA software landing page
  URL: https://ctg.cncr.nl/software/magma
  Why included: official distribution and documentation entry point for the executable and annotations.

## GWAS materials

- GWAS Catalog summary statistics downloads
  URL: https://www.ebi.ac.uk/gwas/downloads/summary-statistics
  Why included: public, curated source for obtaining summary statistics plus study metadata.

- GWAS Catalog summary statistics API docs
  URL: https://www.ebi.ac.uk/gwas/summary-statistics/docs/
  Why included: official programmatic access surface for study and association discovery when the agent needs to fetch data reproducibly.

- GWAS Catalog methods docs
  URL: https://www.ebi.ac.uk/gwas/labs/docs/methods
  Why included: official processing context and curation assumptions for Catalog summary statistics.

## Optional comparator

- scPagwas GitHub repository
  URL: https://github.com/sulab-wmu/scPagwas
  Why included: comparator workflow with a different modeling strategy and different software stack.

- scPagwas paper
  URL: https://doi.org/10.1016/j.xgen.2023.100383
  Why included: method reference for framing optional comparisons to pathway-based GWAS-to-single-cell analysis.

## How to use this inventory

- Prefer the official scDRS docs and repository before third-party tutorials.
- Treat MAGMA gene-level p-values, 10-kb SNP-to-gene mapping, MHC exclusion, and top-1,000-gene selection as the default scDRS runbook unless the task has a documented reason to diverge.
- Summarize downstream results with Benjamini-Hochberg correction, using `FDR < 0.05` for cell-type associations and `FDR < 0.1` for individual significant cells.
- Use GWAS Catalog materials when the task needs public summary statistics or metadata sanity checks.
- Use scPagwas only as a comparator or follow-up, not as a substitute for the requested scDRS workflow.
