# Assets

This directory belongs to the experiment skill `cell-deconvolution`.

The current hardening pass uses deterministic toy or pseudo-data only. No real biological claims should be drawn from these assets.

## What the toy assets validate

- The runnable contract between `examples/`, `scripts/`, `metadata.yaml`, QC sidecars, and the declared deliverables.
- Deterministic computation of gene intersection, NNLS mapping weights, segmentation-aware projection, and held-out prediction QC on tiny synthetic data.
- Repeatable output structure that a general coding agent can regenerate locally without network access.

## What the toy assets do not validate

- Tangram's full constrained optimization, image-derived segmentation, or GPU-backed training loops.
- Method-specific biological correctness on real spatial or single-cell data.
- External database reachability, atlas selection quality, or tissue-specific calibration.

## Upgrade path

- Replace the toy centroids and spot counts with a pinned real benchmark slice that preserves the same deliverable names.
- Swap the NNLS surrogate for Tangram's public `map_cells_to_space` and projection utilities when the heavier dependency stack is acceptable.
- Record atlas choice, segmentation provenance, and model parameters in the report deliverables and `refs.md`.
