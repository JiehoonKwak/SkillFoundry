# Assets

This directory belongs to the experiment skill `squidpy-analysis`.

The current starter uses deterministic toy data only. No real biological claims should be drawn from these assets.

## What the starter validates

- The runnable contract between `examples/`, `scripts/`, `metadata.yaml`, and the declared deliverables.
- A small local path for image-aware spatial neighbors, graph-smoothed label assignment, neighborhood enrichment, co-occurrence, and Moran's I.
- Machine-readable QC files that show how the deliverables were derived from raw toy inputs.

## What the starter approximates

- Squidpy image handling is reduced to precomputed per-spot image summary features instead of an actual `ImageContainer`.
- Neighborhood enrichment is reduced to a deterministic observed-versus-expected score instead of Squidpy's permutation-based z-score.
- Co-occurrence is reduced to a compact near-distance summary rather than the full distance-bin output surface.

## What the starter does not validate

- Statistical significance, multiple-testing correction, or robustness on real spatial omics data.
- Real image cropping, segmentation, or morphology feature extraction from microscopy files.
- Biological interpretation beyond the synthetic contract.

## Upgrade path

- Replace the toy image summaries with pinned public image features or a real Squidpy `ImageContainer` pipeline.
- Carry the same deliverable names and schemas into a richer AnnData-backed implementation.
- Record exact graph, feature, and permutation parameters in `squidpy_report.md` and `refs.md` when promoting the starter.
