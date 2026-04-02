---
name: spatial-deconvolution-portable-skill
description: Use this experiment skill to estimate cell-type abundances in spot-based spatial transcriptomics with a transparent starter choice between cell2location-style, DestVI-style, and baseline models using public, shell-capable file contracts.
---
## Purpose

Estimate cell-type abundances in spot-based spatial transcriptomics with a transparent starter choice between a mapping-style path, a generative-style path, and a simple baseline.

## Source adaptation

This experiment skill is adapted from `source_reference.md`, but it removes private Spatial Agent tools, hidden wrappers, and notebook-only orchestration. The starter stays local, deterministic, and shell-capable.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public Python packages, explicit file contracts, and local synthetic inputs.
- Keep the tiny starter reproducible without network access or GPU requirements.

## Focus terms

`spatial deconvolution, DestVI, cell2location, NNLS`

## Default methods

- `cell2location-like NNLS`
- `DestVI-like smoothed log-NNLS`
- `uniform baseline comparison`

## Recommended resource stack

- `scvi-tools DestVI docs and tutorial`
- `cell2location docs and GitHub`
- `DestVI paper`
- `cell2location paper`
- `spatial deconvolution benchmark or review paper`

## Core workflow

1. Intersect genes between the spatial counts and reference signatures.
2. Build a tiny reference signature matrix from the synthetic cell-type profiles.
3. Run a `cell2location`-shaped deterministic path: library-size normalization plus NNLS abundance fitting.
4. Run a `DestVI`-shaped deterministic path: toy neighbor smoothing, log-scale stabilization, then NNLS abundance fitting.
5. Compare both against a simple baseline using abundance RMSE and correlation, then write the selected abundance matrix and QC artifacts.

## Deliverables

- `cell_type_abundance.h5ad`
- `model_comparison.tsv`
- `deconvolution_summary.md`

## Package scaffold

- `metadata.yaml` defines the experiment contract and required deliverables.
- `examples/toy_input.json` contains only raw synthetic inputs, coordinates, truth tables, and expected invariants.
- `scripts/run_spatial_deconvolution.py` computes the starter outputs from the toy inputs.
- `scripts/run_exercise.py` runs the starter end-to-end and validates the outputs.
- `tests/test_contract.py` and `tests/test_spatial_deconvolution.py` assert structural and numeric invariants.
- `assets/README.md` records what the starter computes versus what remains a surrogate for the full public methods.

## Toy validation

- Run `python3 scripts/run_spatial_deconvolution.py --input examples/toy_input.json --outdir scratch/spatial-deconvolution_toy_run`.
- Run `python3 scripts/run_exercise.py --outdir scratch/spatial-deconvolution_toy_run`.
- Run `python3 -m unittest tests/test_contract.py tests/test_spatial_deconvolution.py`.
- Treat the results as a starter-path contract and QC exercise, not as biological evidence.

## Starter scope

The starter truly computes gene intersection, reference signature assembly, per-spot NNLS abundances, a tiny spatial neighbor graph, a smoothed `DestVI`-style surrogate path, and baseline-vs-model QC metrics on synthetic data. It approximates `cell2location` and `DestVI` with deterministic NNLS plus simple normalization and smoothing instead of probabilistic training, posterior uncertainty, latent state learning, or GPU-backed inference. Atlas selection, tissue priors, uncertainty calibration, and full public-method behavior remain outside scope and should be treated as surrogates in this package.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Record the exact reference atlas, preprocessing decisions, and model parameters in the summary deliverable.
- Stop early when gene overlap, coordinates, or cell-type labels are insufficient for a trustworthy run.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- `model_comparison.tsv` includes the required metric columns from `metadata.yaml`.
- `deconvolution_summary.md` contains the required named sections.
- `cell_type_abundance.h5ad` opens successfully, keeps spatial coordinates, and stores per-spot abundances that sum near 1.
- QC artifacts capture intermediate computations rather than only final files.

## Failure handling

- Stop early and report blockers when feature overlap, positive counts, or coordinates are missing.
- Preserve uncertainty in mixed spots rather than forcing a single confident label.
- Treat the starter outputs as portability and method-shape validation only.

## Hand-off contract

Before finishing, leave the declared deliverables, machine-readable QC files, and a concise run summary so another shell-capable agent can rerun the exact starter path locally.
