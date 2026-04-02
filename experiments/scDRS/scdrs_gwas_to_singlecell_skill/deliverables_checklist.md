# Deliverables Checklist

This checklist defines the exact runtime artifacts and validation gates for a real scDRS execution using the experiment skill.

## Required files

- `01_trait_and_dataset_selection/run_config.yaml`
- `01_trait_and_dataset_selection/selection_report.md`
- `02_gwas_preparation/gwas_input_profile.tsv`
- `02_gwas_preparation/gwas_preflight_summary.json`
- `02_gwas_preparation/trait.sumstats.magma.tsv.gz`
- `03_magma/trait.genes.out`
- `03_magma/magma_run.log`
- `04_scdrs_geneset/trait.magma_gene.zscore.tsv` or an equivalent p-value table
- `04_scdrs_geneset/trait.gs`
- `04_scdrs_geneset/geneset_summary.tsv`
- `05_single_cell_preparation/adata.prepared.h5ad`
- `05_single_cell_preparation/covariates.cov`
- `05_single_cell_preparation/obs_schema.tsv`
- `05_single_cell_preparation/preparation_report.md`
- `06_scdrs_compute_score/trait.score.gz`
- `06_scdrs_compute_score/trait.full_score.gz`
- `06_scdrs_compute_score/compute_score.log`
- `06_scdrs_compute_score/score_summary.tsv`
- `07_scdrs_downstream/downstream_manifest.tsv`
- `07_scdrs_downstream/downstream_report.md`
- `08_positive_control/positive_control_report.md`
- `09_final_application/final_report.md`

## Validation checks

- Dependency check log exists and records whether `python`, `magma`, and `scdrs` were available.
- `run_config.yaml` records trait, dataset, species, source URLs, and intended downstream group columns.
- `gwas_preflight_summary.json` records the original headers, normalized headers, row counts, and blocking issues.
- `gwas_preflight_summary.json` or a linked note records whether the HLA/MHC region was excluded before MAGMA.
- `trait.sumstats.magma.tsv.gz` contains a p-value column and usable variant identifiers.
- `trait.genes.out` is non-empty and MAGMA warnings are preserved in `magma_run.log`.
- `magma_run.log` records that the run used MAGMA gene-level p-values from GWAS summary statistics and whether the annotation resource encoded the default `10 kb` SNP-to-gene window.
- If the annotation resource was rebuilt locally, the run artifacts preserve the corresponding `magma --annotate window=10,10 ...` command or an explicit reason for diverging.
- The MAGMA-to-scDRS reformatted gene-stat table exists and has `GENE` as the first column.
- `trait.gs` is non-empty and `geneset_summary.tsv` records gene-set size.
- `geneset_summary.tsv` records whether the default top `1,000` genes were used or explains the override.
- `covariates.cov` row count matches the number of cells in `adata.prepared.h5ad`.
- `obs_schema.tsv` identifies which `adata.obs` fields were used for grouping or correlation.
- Both `trait.score.gz` and `trait.full_score.gz` exist.
- `score_summary.tsv` reports gene overlap, cell count, `n_ctrl`, and the final `--flag-filter-data` plus `--flag-raw-count` settings.
- `downstream_report.md` explicitly states which downstream analyses ran, which were skipped, and why.
- `downstream_report.md` distinguishes discrete `group-analysis` annotations from continuous numeric `corr-analysis` annotations.
- `downstream_report.md` flags cell groups with fewer than `100` cells as exploratory or low-power, with `<100` cells treated as a practical caution threshold.
- `downstream_report.md` records Benjamini-Hochberg correction and interprets cell-type associations at `FDR < 0.05`.
- `positive_control_report.md` assigns a pass, weak-pass, or fail decision before interpreting the final trait.
- `final_report.md` includes provenance, commands, main findings, unresolved caveats, and the `FDR < 0.1` rule for individual significant cells.

## Partial-run rules

- Missing MAGMA means the run is blocked, not partial.
- Missing scDRS means the run is blocked until installation succeeds.
- Missing biologically meaningful `adata.obs` grouping columns allows a partial run only for `compute-score`; downstream grouping must be marked incomplete.
- Failed positive control means the final trait interpretation is exploratory only and must be labeled unvalidated.
