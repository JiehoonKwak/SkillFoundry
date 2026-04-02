# scDRS Execution Checklist

Fill this during a real run. Replace placeholders rather than deleting fields.

## Run identity

- Run date:
- Operator or agent:
- Trait name:
- Trait slug:
- GWAS source URL or accession:
- Single-cell dataset path:
- Species:

## Dependency check

- `python` available:
- `scdrs` available:
- `magma` available:
- Dependency log path:
- If blocked, exact blocker:

## Stage 1: Trait and dataset selection

- Biological rationale:
- Expected enriched cells or states:
- Positive-control trait selected:
- `run_config.yaml` written:
- `selection_report.md` written:

## Stage 2: GWAS preparation

- Original GWAS filename:
- Original headers:
- Normalized headers:
- Genome build assumption:
- Ancestry assumption:
- Sample-size field used:
- HLA/MHC region exclusion applied:
- Blocking format issues:
- `gwas_preflight_summary.json` written:
- `trait.sumstats.magma.tsv.gz` written:

## Stage 3: MAGMA

- LD reference used:
- Gene annotation used:
- Annotation command recorded if rebuilt locally:
- SNP-to-gene window used:
- Default `10 kb` window kept:
- MAGMA command recorded:
- `trait.genes.out` written:
- Top 5 genes checked:
- MAGMA warnings:
- If failed, blocker note path:

## Stage 4: `scdrs munge-gs`

- Final command:
- Reformatted MAGMA table path:
- Statistic mode: z-score / p-value
- Input statistic column used:
- Default top `1,000` genes kept:
- Output gene-set size:
- `trait.gs` written:
- `geneset_summary.tsv` written:

## Stage 5: Single-cell preparation

- Input `.h5ad` shape:
- `adata.obs` columns reviewed:
- Count-state assumption:
- Covariates used:
- Group columns available:
- Correlation columns available:
- `adata.prepared.h5ad` written:
- `covariates.cov` written:
- `obs_schema.tsv` written:

## Stage 6: `scdrs compute-score`

- Final command:
- `--flag-filter-data` value:
- `--flag-raw-count` value:
- Why those flags match the prepared `h5ad`:
- `--n-ctrl` value:
- `trait.score.gz` written:
- `trait.full_score.gz` written:
- Gene overlap count:
- Main warnings:

## Stage 7: `scdrs perform-downstream`

- Group analyses run:
- Correlation analyses run:
- continuous numeric columns used for `corr-analysis`:
- Small groups flagged as exploratory (`<100` cells):
- Gene analysis run:
- BH correction recorded:
- Cell-type significance rule (`FDR < 0.05`) applied:
- Skipped analyses and reasons:
- `downstream_manifest.tsv` written:
- `downstream_report.md` written:

## Stage 8: Positive control

- Positive-control trait:
- Expected biology:
- Observed biology:
- Decision: pass / weak-pass / fail
- Follow-up if not pass:
- `positive_control_report.md` written:

## Stage 9: Final disease-trait application

- Main enriched groups:
- Main cell-state correlations:
- Individual significant cell rule (`FDR < 0.1`) applied:
- Main caveats:
- Validation status:
- `final_report.md` written:

## Sign-off

- All required files present:
- Any blocked stages:
- Ready for review:
