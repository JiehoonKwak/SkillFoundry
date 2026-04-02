---
name: gene-imputation-portable-skill
description: Use this experiment skill to Impute unmeasured or sparse genes in spatial data from a paired scRNA-seq reference using portable workflows and explicit validation with public tools and explicit file contracts that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Impute unmeasured or sparse genes in spatial data from a paired scRNA-seq reference using portable workflows and explicit validation.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`gene imputation, spatial, tangram`

## Default methods

- Tangram projection
- held-out gene validation
- marker recovery checks

## Recommended resource stack

- Tangram documentation and tutorials
- Tangram GitHub repository
- Tangram Nature Methods paper
- Squidpy and AnnData documentation for portable I/O and spatial containers

## Core workflow

1. Confirm reference-to-spatial feature compatibility before mapping or projection.
2. Aggregate the toy scRNA-seq reference into cell-type centroids over the shared and held-out genes.
3. Estimate spot-to-centroid weights with simplex-constrained NNLS and save cosine-similarity QC alongside the fitted weights.
4. Project held-out genes from the fitted weights, score marker recovery, and record assignment entropy before writing the final deliverables.

## Deliverables

- `imputed_expression.h5ad`
- `heldout_metrics.tsv`
- `imputation_report.md`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract for `gene-imputation`.
- `examples/toy_input.json` contains raw synthetic reference cells, spatial shared genes, held-out truth, marker sets, and expected invariants.
- `scripts/run_gene_imputation.py` computes centroid mapping weights, held-out predictions, and machine-readable QC without network access.
- `scripts/validate_outputs.py` rechecks the declared deliverables plus task-specific numeric invariants.
- `tests/test_gene_imputation.py` asserts mapping-weight, metric, and AnnData structure invariants.
- `assets/README.md` explains the starter boundary and what remains outside scope.

## Toy validation

- Run `python3 scripts/run_gene_imputation.py --input examples/toy_input.json --outdir scratch/gene-imputation_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/gene-imputation_toy_run --input examples/toy_input.json` to recheck deliverables and toy invariants.
- Run `python3 -m unittest tests/test_contract.py tests/test_gene_imputation.py` for the skill-local test pass.
- Inspect `intermediate_qc.json` for cosine similarities, mapping weights, spot assignments, and marker recovery.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact atlas, reference slice, model version, normalization choices, and solver settings in the final report or provenance notes.
- Keep stop conditions explicit when feature overlap, spatial metadata, or reference labels are insufficient for a trustworthy run.

## Starter scope

- Computed: gene intersection, cell-type centroids, simplex-constrained NNLS mapping weights, projected spot-by-gene expression, held-out metrics, marker recovery, and assignment entropy on a tiny synthetic problem.
- Approximated: full Tangram optimization is replaced by centroid-level deterministic fitting rather than PyTorch cell-level training.
- Surrogate: held-out truth is synthetic and bundled with the toy input, so this starter validates the workflow shape and QC logic rather than real biological performance.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- `heldout_metrics.tsv` keeps the required columns in `metadata.yaml`.
- `imputation_report.md` contains the named sections used for downstream review.
- `imputed_expression.h5ad` opens successfully, contains mapping weights, and stores non-empty projected expression.
- `intermediate_qc.json` captures the machine-readable intermediate checks used by the tests.

## Failure handling

- Stop early and document blockers when core metadata, coordinates, identifiers, or feature overlap are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the intermediate QC JSON, a concise run summary, and enough provenance for another shell-capable agent to rerun the experiment.
