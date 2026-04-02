# References

- scDRS documentation home
  URL: https://martinjzhang.github.io/scDRS/
  Note: official package overview, install notes, and usage entry points.

- scDRS CLI reference
  URL: https://martinjzhang.github.io/scDRS/reference_cli.html
  Note: canonical command surface for `munge-gs`, `compute-score`, and `perform-downstream`.

- scDRS file formats
  URL: https://martinjzhang.github.io/scDRS/file_format.html
  Note: central file contracts for `.gs`, `.cov`, p-value and z-score inputs, and downstream output naming.

- scDRS FAQ
  URL: https://martinjzhang.github.io/scDRS/faq.html
  Note: concise statement of the standard end-to-end workflow, the top-1,000-gene convention, and downstream-analysis expectations.

- scDRS GitHub repository
  URL: https://github.com/martinjzhang/scDRS
  Note: code repository for version inspection, examples, tests, and issue tracking when CLI behavior differs by release.

- scDRS issue #2: custom gene set
  URL: https://github.com/martinjzhang/scDRS/issues/2
  Note: maintainer-linked discussion quoting the paper's standard MAGMA setup: gene-level p-values from GWAS summary statistics, 10-kb SNP-to-gene window, and top 1,000 genes.

- Zhang, Hou et al. 2022. Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC9891382/
  Note: method paper defining scDRS, the MHC exclusion, the standard MAGMA-to-top-1,000-gene workflow, and the reported FDR interpretation thresholds.

- Zhang, Hou et al. 2022. Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data
  URL: https://www.nature.com/articles/s41588-022-01167-z
  Note: publisher version of the primary method paper.

- de Leeuw et al. 2015. MAGMA: Generalized Gene-Set Analysis of GWAS Data
  URL: https://doi.org/10.1371/journal.pcbi.1004219
  Note: canonical MAGMA method reference for SNP-to-gene aggregation and LD-aware gene analysis.

- MAGMA software landing page
  URL: https://ctg.cncr.nl/software/magma
  Note: official distribution/documentation entry point for binaries, annotations, and manuals.

- GWAS Catalog summary statistics downloads
  URL: https://www.ebi.ac.uk/gwas/downloads/summary-statistics
  Note: canonical public source for published and pre-published GWAS summary statistics plus metadata.

- GWAS Catalog summary statistics API documentation
  URL: https://www.ebi.ac.uk/gwas/summary-statistics/docs/
  Note: official API and field-surface reference for programmatic summary-statistics discovery and retrieval.

- GWAS Catalog methods overview
  URL: https://www.ebi.ac.uk/gwas/labs/docs/methods
  Note: official documentation hub that links curation and summary-statistics processing details.

- scPagwas GitHub repository
  URL: https://github.com/sulab-wmu/scPagwas
  Note: optional comparator workflow based on pathway-activity transformation in R/Seurat rather than scDRS.

- Ma et al. 2023. Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data
  URL: https://doi.org/10.1016/j.xgen.2023.100383
  Note: scPagwas method paper for optional comparison framing.
