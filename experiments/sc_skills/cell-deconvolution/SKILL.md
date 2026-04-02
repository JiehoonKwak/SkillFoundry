---
name: cell-deconvolution-portable-skill
description: Use this experiment skill to Deconvolve spot-based spatial measurements to likely single-cell assignments using mapping constraints and explicit QC with public tools and explicit file contracts that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Deconvolve spot-based spatial measurements to likely single-cell assignments using mapping constraints and explicit QC.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`spatial, deconvolution, tangram`

## Default methods

- Tangram constrained mapping
- segmentation-aware projection
- mapping QC

## Recommended resource stack

- `Tangram documentation`
- `Tangram GitHub repository`
- `Squidpy documentation`
- `Scanpy documentation`
- `Tangram paper`
- `Spatial integration benchmark`
- `CELLxGENE Census documentation`
- `CellTypist documentation`
- `Azimuth documentation`
- `SpatialData documentation`
- `CellTypist documentation`

## Core workflow

1. Intersect the reference and spatial genes, then split them into fit genes and held-out QC genes.
2. Estimate Tangram-style spot-to-reference weights with deterministic NNLS over fit-gene centroids.
3. Project those weights onto per-spot segmentation counts to approximate constrained single-cell assignments.
4. Save mapping QC, including held-out gene prediction, assignment entropy, and marker recovery.

## Deliverables

- `mapping_matrix.h5ad`
- `cell_assignment.tsv`
- `deconvolution_report.md`

## Package scaffold

- `metadata.yaml` defines the toy contract and required deliverables for `cell-deconvolution`.
- `examples/toy_input.json` provides raw synthetic centroids, spot counts, coordinates, segmentation counts, and expected invariants.
- `scripts/run_cell_deconvolution.py` computes NNLS weights, segment projection, held-out QC, and the declared deliverables.
- `tests/test_cell_deconvolution.py` checks numeric and structural invariants of the runnable starter.
- `assets/README.md` explains what the pseudo-data validates and what it does not.

## Toy validation

- Run `python3 scripts/run_cell_deconvolution.py --input examples/toy_input.json --outdir scratch/cell-deconvolution_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 -m unittest discover -s tests` for the skill-local starter checks.
- Treat these outputs as toy or pseudo-data validation of the file contract, not as biological evidence.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable interface.
- Pin the exact atlas, database snapshot, model version, and key CLI parameters in the final report or provenance notes.
- Keep stop conditions explicit when feature overlap, spatial metadata, or reference labels are insufficient for a trustworthy run.

## Starter scope

- Computed: gene intersection, log1p normalization, NNLS mapping weights, held-out gene prediction, assignment entropy, and integerized segment-count projection.
- Approximated: Tangram constrained mode is replaced with centroid-level NNLS plus provided segmentation counts instead of learned cell-to-space optimization.
- Surrogate only: image-derived segmentation, GPU training, atlas choice, and real-tissue biological interpretation remain outside starter scope.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- TSV outputs keep the required columns in `metadata.yaml`.
- Markdown reports contain the named sections used for downstream review.
- `.h5ad` and pseudo-`.h5mu` outputs open successfully and are non-empty.

## Failure handling

- Stop early and document blockers when core metadata, coordinates, identifiers, or feature overlap are missing.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve unresolved ambiguity rather than inventing certainty.

## Hand-off contract

Before finishing, leave behind the declared deliverables, a concise run summary, and enough provenance for another shell-capable agent to rerun the experiment.
