---
name: scdrs-gwas-to-singlecell-skill
description: Use this skill to run a portable scDRS workflow that maps GWAS disease-risk signals onto single-cell RNA-seq data using explicit MAGMA, scDRS, AnnData, and covariate file contracts without private tool dependencies.
---

## Purpose

Provide a strict, reusable workflow for taking a disease or trait GWAS from summary statistics to single-cell disease-relevance scoring with scDRS. The workflow is written for general coding agents and makes every dependency, file contract, checkpoint, and failure mode explicit.

## When to use

- You need a disease-risk enrichment workflow on scRNA-seq data with a public, scriptable CLI path.
- You want a stepwise plan for MAGMA plus scDRS that works outside Spatial Agent.
- You need reproducible output artifacts that can be reviewed, rerun, or handed to another agent.

## When not to use

- You only need GWAS QC or only need single-cell preprocessing.
- You do not have access to legal, appropriately licensed GWAS summary statistics or single-cell data.
- Your single-cell input is not RNA-based and you are not prepared to adapt the workflow with custom controls.

## Required capabilities

- Shell access with standard Unix utilities.
- Python environment that can run `python`, `pip`, and read/write `h5ad`.
- A working MAGMA installation or permission to install/download it.
- A working scDRS installation or permission to install it.

## Dependency check

Run this before doing any data work:

```bash
mkdir -p work/scdrs_run/logs
{
  echo "python: $(command -v python || true)"
  echo "pip: $(command -v pip || true)"
  echo "magma: $(command -v magma || true)"
  echo "scdrs: $(command -v scdrs || true)"
  python --version 2>&1 || true
  pip show scdrs 2>/dev/null || true
} | tee work/scdrs_run/logs/dependency_check.txt
```

If `scdrs` is missing:

```bash
pip install scdrs
python -c "import scdrs; print(scdrs.__version__)"
```

If `magma` is missing:

1. Stop and record the blocker in the run log.
2. Obtain MAGMA from the official documentation or licensed institutional distribution.
3. Add the binary to `PATH` and rerun `magma --help`.

Do not pretend MAGMA ran if it is absent. The workflow cannot produce valid scDRS GWAS-derived gene statistics without it.

## Workspace contract

Use a run directory with stable stage folders:

```text
work/scdrs_run/
  inputs/
    gwas/
    single_cell/
    refs/
  01_trait_and_dataset_selection/
  02_gwas_preparation/
  03_magma/
  04_scdrs_geneset/
  05_single_cell_preparation/
  06_scdrs_compute_score/
  07_scdrs_downstream/
  08_positive_control/
  09_final_application/
  logs/
```

Keep every intermediate file. Do not overwrite silently.

## Input contract

### GWAS input

Provide a GWAS summary-statistics file under `inputs/gwas/` plus a short metadata note containing:

- trait name and slug
- source URL or accession
- genome build if known
- ancestry if known
- sample size or effective sample size
- effect columns present
- p-value column present
- whether the file is raw or harmonized

Minimum fields expected for MAGMA preparation:

- SNP identifier or `CHR` + `BP`
- p-value
- sample size or effective sample size
- effect allele and non-effect allele if you plan to normalize/harmonize before MAGMA

Preferred additional fields:

- beta or odds ratio
- standard error
- INFO
- EAF

If the GWAS format is incompatible:

1. Write `02_gwas_preparation/gwas_format_issue.md` describing missing columns and exact headers found.
2. Attempt a deterministic header-normalization pass.
3. If core identifiers or p-values are missing after normalization, stop and request a compatible GWAS file.

### Single-cell input

Provide one `.h5ad` under `inputs/single_cell/` with:

- cells in `adata.obs_names`
- genes in `adata.var_names`
- expression matrix in `adata.X` or `.raw`
- enough metadata in `adata.obs` to support downstream grouping or covariate design

Preferred `adata.obs` fields:

- sample or donor ID
- batch
- cell type
- tissue or region
- disease or condition label if relevant
- QC metrics if available

If `adata.obs` lacks needed annotations:

1. Save `05_single_cell_preparation/obs_schema.tsv` with all observed columns and dtypes.
2. Add only portable, derivable annotations first, such as `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, or a copied batch/sample column.
3. If no biologically meaningful grouping exists for downstream analysis, still run `compute-score` but mark `perform-downstream` as partial and explain the limitation in `07_scdrs_downstream/downstream_blockers.md`.

## Output contract

Always materialize these stable files, even if the underlying tool emits differently named files:

- `01_trait_and_dataset_selection/run_config.yaml`
- `01_trait_and_dataset_selection/selection_report.md`
- `02_gwas_preparation/gwas_input_profile.tsv`
- `02_gwas_preparation/gwas_preflight_summary.json`
- `02_gwas_preparation/trait.sumstats.magma.tsv.gz`
- `03_magma/trait.genes.out`
- `03_magma/magma_run.log`
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

## Stage 1: Trait and dataset selection

Pick a trait-dataset pair that is biologically plausible and statistically powered.

Selection rules:

- Prefer GWAS with heritability z-score greater than 5 or total sample size greater than 100K when available.
- Prefer single-cell datasets whose tissue, disease context, or dominant cell compartments plausibly match the GWAS trait.
- Record why the pair is a fit and any known caveats.

Write `run_config.yaml` with at least:

```yaml
trait_name: "example_trait"
trait_slug: "example_trait"
gwas_source_url: "https://..."
gwas_build: "GRCh37"
gwas_ancestry: "EUR"
single_cell_h5ad: "work/scdrs_run/inputs/single_cell/example.h5ad"
h5ad_species: "human"
gs_species: "human"
positive_control_trait: "autoimmune_trait_or_known_reference"
group_columns:
  - cell_type
  - batch
cell_score_columns: []
gene_score: true
```

Write `selection_report.md` summarizing:

- why this GWAS was selected
- why this dataset was selected
- expected positive-control biology
- known power or compatibility risks

## Stage 2: GWAS summary-statistics preparation

Profile and normalize the GWAS before MAGMA.

scDRS's standard workflow starts from MAGMA gene-level p-values computed from GWAS summary statistics. The paper and scDRS issue discussion describe a default setup that maps SNPs to genes with a 10-kb window around the gene body and then selects the top 1,000 genes by MAGMA p-value as putative disease genes. Keep these defaults explicit in the run log even when you intentionally change them.

Checklist:

1. Confirm delimiter, compression, row count, and header names.
2. Normalize common aliases into a stable table with columns such as `SNP`, `CHR`, `BP`, `A1`, `A2`, `P`, `N`, `BETA`, `OR`, `SE`.
3. Exclude the HLA/MHC region before MAGMA or document why you cannot. The scDRS paper removed the MHC region because of its unusual LD and genetic architecture, so this should be the default behavior rather than an afterthought.
4. Drop rows with invalid p-values, impossible chromosomes, or missing identifiers.
5. Record allele and build assumptions.
6. Save a machine-readable summary.

Example profiling commands:

```bash
zcat -f work/scdrs_run/inputs/gwas/trait.sumstats.tsv.gz | head -5
zcat -f work/scdrs_run/inputs/gwas/trait.sumstats.tsv.gz | wc -l
```

Checkpoint:

- `gwas_preflight_summary.json` records headers, row counts, missingness, and blocking issues.
- `trait.sumstats.magma.tsv.gz` has the exact columns you will feed to MAGMA.

Stop if:

- no valid p-value column remains
- no SNP ID or no chromosome-position combination remains
- the file mixes genome builds and you cannot repair it deterministically

## Stage 3: MAGMA gene-level statistics

Run MAGMA only after the GWAS file is normalized and the required reference inputs are present.

Required reference inputs under `inputs/refs/`:

- LD reference panel compatible with the GWAS ancestry
- MAGMA gene-location annotation resources for the target genome build

Recommended default MAGMA configuration:

- start from the normalized GWAS p-value table in `02_gwas_preparation/trait.sumstats.magma.tsv.gz`
- compute MAGMA gene-level p-values from summary statistics
- use a gene annotation file that was built with a 10-kb SNP-to-gene mapping window around the gene body, or explicitly document a different window
- record whether the 10-kb window was baked into the annotation resource you used

If you need to build the MAGMA annotation resource yourself, do that before the gene-analysis command:

```bash
magma --annotate \
  window=10,10 \
  --snp-loc work/scdrs_run/inputs/refs/ld_reference_prefix.bim \
  --gene-loc work/scdrs_run/inputs/refs/gene_locations.tsv \
  --out work/scdrs_run/inputs/refs/gene_annotation
```

MAGMA accepts a PLINK `.bim` file as the SNP-location input for `--snp-loc`. If you use a different text file, make sure it follows MAGMA's SNP-location contract.

Then feed the resulting `work/scdrs_run/inputs/refs/gene_annotation.genes.annot` file into the gene-analysis step below.

Example MAGMA pattern:

```bash
magma \
  --bfile work/scdrs_run/inputs/refs/ld_reference_prefix \
  --pval work/scdrs_run/02_gwas_preparation/trait.sumstats.magma.tsv.gz use=SNP,P ncol=N \
  --gene-model snp-wise=mean \
  --gene-annot work/scdrs_run/inputs/refs/gene_annotation.genes.annot \
  --out work/scdrs_run/03_magma/trait
```

If your GWAS lacks `N`, replace `ncol=N` with the correct MAGMA argument for the available sample-size column or stop and document the mismatch.

Checkpoint:

- `trait.genes.out` exists and is non-empty.
- `magma_run.log` captures the command, MAGMA version, reference panel, and any warnings.
- The top gene-level hits are plausible and not all missing or identical.
- The log states whether the MHC exclusion was applied and whether the gene annotation resource already encoded the default 10-kb mapping window.

If MAGMA fails:

1. Save stderr/stdout to `magma_run.log`.
2. Check the reference panel prefix, annotation file path, genome build, and `use=` column mapping.
3. If the issue is unresolved after one corrected retry, stop and write `03_magma/magma_blocker.md`.

## Stage 4: Convert MAGMA output to scDRS gene sets

Use `scdrs munge-gs` only after reformatting the MAGMA output into an scDRS-friendly gene-by-trait table whose first column is `GENE` and whose remaining trait columns are numeric gene-level statistics. For the standard scDRS workflow, prefer MAGMA gene-level p-values and select the top 1,000 genes as putative disease genes unless you have a documented reason to change `--n-max`.

Preferred pattern when starting from MAGMA gene-level p-values:

```bash
scdrs munge-gs \
  --pval-file work/scdrs_run/04_scdrs_geneset/trait.magma_gene.pval.tsv \
  --out-file work/scdrs_run/04_scdrs_geneset/trait.gs \
  --n-max 1000
```

Alternative pattern when MAGMA z-scores are available:

```bash
scdrs munge-gs \
  --zscore-file work/scdrs_run/04_scdrs_geneset/trait.magma_gene.zscore.tsv \
  --out-file work/scdrs_run/04_scdrs_geneset/trait.gs \
  --n-max 1000 \
  --weight zscore
```

If you only have gene-level p-values after MAGMA post-processing, switch to the p-value input mode for your installed scDRS release and keep the same output contract.

Before running the command:

1. Reformat `trait.genes.out` into a clean TSV with `GENE` as the first column.
2. Use one trait statistic column per trait.
3. Prefer the MAGMA p-value column when available; use z-scores only when that is the documented input path for your scDRS release.
4. Keep `--n-max 1000` as the default and record any deviation explicitly.
5. Record the exact command from `scdrs munge-gs --help` in the log because option names can drift across releases.

Checkpoint:

- `trait.gs` exists and is non-empty.
- `trait.magma_gene.zscore.tsv` or the equivalent p-value table exists and is archived.
- `geneset_summary.tsv` lists trait name, gene count, weighting mode, and any genes dropped.
- Gene-set size is moderate. Very small sets or sets covering a large fraction of all genes should be flagged as weak or invalid.

If `munge-gs` cannot parse the MAGMA output:

1. Inspect the columns in `trait.genes.out`.
2. Reformat to the scDRS-required p-value or z-score file contract.
3. If you still cannot identify a gene-name column plus a numeric statistic column, stop and document the incompatibility.

## Stage 5: Prepare `.h5ad` and `.cov`

Inspect the AnnData object before running scDRS.

Example inspection:

```bash
python - <<'PY'
import scanpy as sc
adata = sc.read_h5ad("work/scdrs_run/inputs/single_cell/example.h5ad")
print(adata)
print(list(adata.obs.columns))
print(list(adata.var.columns)[:20])
print(adata.obs.head())
PY
```

Requirements:

- `adata.obs_names` must uniquely identify cells.
- Gene names must be meaningful for the chosen species and compatible with the gene set.
- Decide whether `adata.X` contains raw counts. If not, set the scDRS flags accordingly.

Create `covariates.cov` as a tab-separated table indexed by `adata.obs_names`.

Recommended covariates:

- batch
- donor or sample
- total counts
- number of detected genes
- mitochondrial percentage

If no external covariates are available, derive QC covariates from the AnnData object and document them.

Checkpoint:

- `adata.prepared.h5ad` is saved after any required normalization or metadata cleanup.
- `covariates.cov` row count matches `adata.n_obs`.
- `obs_schema.tsv` lists every `adata.obs` column, dtype, missingness, and whether it is suitable for grouping or correlation tests.

## Stage 6: Run `scdrs compute-score`

The CLI documentation uses dash-style flags such as `--flag-filter-data`; scDRS treats dash and underscore spellings as equivalent, but prefer the documented dash form in this runbook to reduce confusion.

Example pattern:

```bash
scdrs compute-score \
  --h5ad-file work/scdrs_run/05_single_cell_preparation/adata.prepared.h5ad \
  --h5ad-species human \
  --cov-file work/scdrs_run/05_single_cell_preparation/covariates.cov \
  --gs-file work/scdrs_run/04_scdrs_geneset/trait.gs \
  --gs-species human \
  --out-folder work/scdrs_run/06_scdrs_compute_score \
  --flag-filter-data False \
  --flag-raw-count False \
  --n-ctrl 1000
```

Set `--flag-raw-count` and `--flag-filter-data` to match the actual `.h5ad` state. `--flag-filter-data` controls whether scDRS applies its own minimal cell and gene filtering. `--flag-raw-count` controls whether scDRS should treat `adata.X` as raw counts and perform size-factor normalization plus `log1p` internally. Do not copy the example values blindly.

Practical rule:

- If `adata.X` still contains raw counts and you want scDRS to handle normalization, use `--flag-raw-count True`.
- If `adata.X` is already normalized or log-transformed, use `--flag-raw-count False`.
- If you already performed the same minimal filtering externally and want to preserve the prepared object as-is, use `--flag-filter-data False`.

Checkpoint:

- `trait.score.gz` exists.
- `trait.full_score.gz` exists.
- `score_summary.tsv` records cell count, gene overlap count, trait name, `n_ctrl`, and the distribution of raw and normalized scores.

If scDRS fails:

1. Check species mismatch between `.h5ad` genes and `.gs`.
2. Check that `.cov` uses the same cell IDs as `adata.obs_names`.
3. Check whether `adata.X` count state matches `--flag-raw-count`.
4. Fall back to the Python API only if the CLI is unavailable but the package imports cleanly; record the fallback explicitly.

## Stage 7: Run `scdrs perform-downstream`

Use the full score file plus meaningful `adata.obs` annotations.

Example pattern:

```bash
scdrs perform-downstream \
  --h5ad-file work/scdrs_run/05_single_cell_preparation/adata.prepared.h5ad \
  --score-file work/scdrs_run/06_scdrs_compute_score/trait.full_score.gz \
  --out-folder work/scdrs_run/07_scdrs_downstream \
  --group-analysis cell_type \
  --corr-analysis total_counts \
  --gene-analysis
```

Use only annotations that actually exist in `adata.obs` and have interpretable cardinality.

Interpretation rules:

- `--group-analysis` is for discrete cell-group annotations such as cell type, tissue, condition, or donor-derived clusters.
- `--corr-analysis` is only for continuous numeric annotations in `adata.obs`, such as pseudotime, QC metrics, module scores, or experimentally measured gradients. Do not pass arbitrary labels, string categories, or encoded IDs.
- Treat group-level results as more reliable when each tested group has adequate size. As a practical heuristic, flag groups with fewer than 100 cells as exploratory or low-power instead of treating them as equally stable.

Checkpoint:

- `downstream_manifest.tsv` lists every downstream result file and the annotation or analysis mode used.
- `downstream_report.md` explains which group, correlation, and gene analyses succeeded and which were skipped.
- `downstream_report.md` records the multiple-testing rule used for summaries. Use Benjamini-Hochberg correction when summarizing downstream result tables, interpret cell-type associations at `FDR < 0.05`, and interpret individual significant cells at `FDR < 0.1`.

If `adata.obs` lacks grouping annotations:

- Skip group analysis.
- Keep correlation analysis only for valid numeric columns.
- Document the missing group labels as a partial-run limitation.

## Stage 8: Positive-control validation

Before claiming success on the target disease trait, run one positive-control analysis with a trait-dataset pair that should produce known biology.

Examples of positive-control logic:

- immune-mediated GWAS on PBMC or immune-rich tissue
- neuropsychiatric GWAS on brain single-cell data
- lipid or liver-associated GWAS on hepatocyte-rich data

Validation rules:

- The chosen positive control must have a clear biological expectation written in advance.
- The report must state whether the expected cell type or state appears among the strongest associations.
- If the positive control fails, do not over-interpret the target-trait result. Diagnose data compatibility, gene overlap, and covariate choices first.

Write `08_positive_control/positive_control_report.md` with:

- positive-control trait and source
- expected enriched cell populations
- observed result summary
- pass, weak-pass, or fail decision
- next action if not pass

## Stage 9: Final disease-trait application

Only after the positive control is acceptable:

1. Run the target trait through the same fixed pipeline.
2. Summarize score distributions, top group-level enrichments, and notable cell-state correlations using Benjamini-Hochberg-adjusted summaries.
3. Record all caveats: power, ancestry mismatch, missing annotations, gene overlap, covariates, and unresolved warnings.

Write `09_final_application/final_report.md` with:

- data provenance
- commands run
- key outputs
- top biological findings
- significance interpretation, including `FDR < 0.05` for cell-type associations and `FDR < 0.1` for individual cells
- validation status
- limitations
- exact rerun paths

## Portable validation logic

A run is minimally acceptable only if all of the following are true:

- dependencies were checked and versions were recorded
- GWAS provenance and compatibility were documented
- MAGMA produced a non-empty gene-level result
- scDRS produced both `.score.gz` and `.full_score.gz`
- downstream analysis was either executed or explicitly marked partial with a reason
- a positive-control decision was recorded before the final interpretation

## Failure-handling summary

- Missing `magma`: install or acquire MAGMA, then stop until `magma --help` works.
- Missing `scdrs`: install with `pip install scdrs`, confirm import, and rerun.
- Incompatible GWAS format: normalize headers once; if core identifiers or p-values remain absent, stop and request a new file.
- Missing `adata.obs` annotations: generate schema and QC-derived covariates, run `compute-score`, and mark downstream grouping as partial if biologically meaningful labels cannot be recovered.
- Low gene overlap between `.gs` and `.h5ad`: check species, identifier type, capitalization, and build assumptions before interpreting results.
- Weak or null positive control: treat the final trait run as unvalidated.

## Provenance

- scDRS docs: https://martinjzhang.github.io/scDRS/
- scDRS CLI docs: https://martinjzhang.github.io/scDRS/reference_cli.html
- scDRS file formats: https://martinjzhang.github.io/scDRS/file_format.html
- scDRS FAQ: https://martinjzhang.github.io/scDRS/faq.html
- scDRS paper: https://www.nature.com/articles/s41588-022-01167-z
- MAGMA paper: https://doi.org/10.1371/journal.pcbi.1004219
- GWAS Catalog summary statistics: https://www.ebi.ac.uk/gwas/downloads/summary-statistics
- scPagwas comparator: https://github.com/sulab-wmu/scPagwas
