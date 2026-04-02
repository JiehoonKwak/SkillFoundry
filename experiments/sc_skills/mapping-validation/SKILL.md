---
name: mapping-validation-portable-skill
description: Use this experiment skill to validate spatial mapping outputs with quantitative metrics, biological sanity checks, and hold-out experiments using portable local workflows that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Validate spatial mapping outputs with quantitative metrics, biological sanity checks, and hold-out experiments.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes private Spatial Agent tooling, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`mapping validation, tangram, spatial`

## Default methods

- cross-validation
- marker recovery
- spatial coherence checks

## Recommended resource stack

- Tangram documentation and tutorial
- Tangram utility docs for `compare_spatial_geneexp`, `eval_metric`, and `cross_val`
- Tangram GitHub repository
- Tangram Nature Methods paper
- Spatial integration benchmark
- Squidpy paper or docs for neighborhood-graph upgrade paths

## Core workflow

1. Confirm shared-gene overlap between the reference and spatial inputs.
2. Aggregate the toy scRNA-seq reference into cell-type centroids and score cosine similarity to each spatial spot.
3. Fit simplex-constrained NNLS mapping weights from centroids to spatial spots.
4. Validate those weights with held-out gene prediction, leave-one-gene-out cross-validation, marker recovery, assignment entropy, and a tiny nearest-neighbor spatial coherence check.

## Deliverables

- `validation_metrics.tsv`
- `holdout_summary.md`
- `mapping_validation_report.md`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract for `mapping-validation`.
- `examples/toy_input.json` contains raw synthetic reference cells, spatial spots, marker sets, coordinates, and expected invariants.
- `scripts/run_mapping_validation.py` computes centroid mapping weights, held-out metrics, leave-one-gene-out validation, and machine-readable QC without network access.
- `scripts/validate_outputs.py` rechecks the declared deliverables plus task-specific numeric invariants.
- `tests/test_mapping_validation.py` asserts mapping-weight, metric, and QC-table invariants.
- `assets/README.md` explains the starter boundary and what remains outside scope.

## Toy validation

- Run `python3 scripts/run_mapping_validation.py --input examples/toy_input.json --outdir scratch/mapping-validation_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/mapping-validation_toy_run --input examples/toy_input.json` to recheck deliverables and toy invariants.
- Run `python3 -m unittest tests/test_contract.py tests/test_mapping_validation.py` for the skill-local test pass.
- Inspect `intermediate_qc.json`, `gene_cv_qc.tsv`, `spot_assignment_qc.tsv`, and `spatial_coherence_qc.tsv` for the computed intermediate checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact atlas slice, normalization choices, model version, and solver settings in the final report or provenance notes.
- Keep stop conditions explicit when feature overlap, spatial metadata, or reference labels are insufficient for a trustworthy run.

## Starter scope

- Computed: gene intersection, centroid-level cosine similarity, simplex-constrained NNLS mapping weights, held-out gene prediction, leave-one-gene-out cross-validation, marker recovery, assignment entropy, and a tiny nearest-neighbor coherence check on synthetic data.
- Approximated: full Tangram cell-level optimization and full Squidpy graph statistics are replaced by deterministic centroid fitting plus a lightweight local graph heuristic.
- Surrogate: held-out truth and biological sanity rules are bundled synthetic targets, so the starter validates workflow shape and QC logic rather than real biological performance.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- `validation_metrics.tsv` keeps the required columns in `metadata.yaml`.
- Markdown deliverables contain the named sections used for downstream review.
- `intermediate_qc.json`, `gene_cv_qc.tsv`, `spot_assignment_qc.tsv`, and `spatial_coherence_qc.tsv` stay consistent with the fitted mapping weights.

## Failure handling

- Stop early and document blockers when core metadata, coordinates, identifiers, or feature overlap are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the intermediate QC files, a concise run summary, and enough provenance for another shell-capable agent to rerun the experiment.
