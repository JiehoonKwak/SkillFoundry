# Assets

This directory belongs to the experiment skill `mapping-validation`.

The current starter uses deterministic synthetic inputs only. No real biological claims should be drawn from these assets.

## Starter boundary

- Computed locally: shared-gene intersection, reference centroids, cosine similarity, simplex-constrained NNLS mapping weights, held-out prediction, leave-one-gene-out cross-validation, marker recovery, assignment entropy, and a tiny nearest-neighbor coherence check.
- Approximated: full Tangram optimization is replaced by centroid-level deterministic fitting, and richer Squidpy neighborhood statistics are replaced by a minimal graph heuristic.
- Surrogate: held-out truth and biological sanity expectations are synthetic, so the starter validates workflow shape and QC logic rather than public benchmark performance.

## What the starter validates

- The runnable contract between `examples/`, `scripts/`, `metadata.yaml`, and the declared deliverables.
- Machine-readable QC artifacts that another shell-capable agent can inspect or test.
- Deterministic numerical behavior for mapping weights, held-out metrics, and graph-based coherence on tiny toy data.

## What the starter does not validate

- Biological correctness on real spatial transcriptomics datasets.
- Tangram training dynamics, PyTorch optimization behavior, or cell-level mapping outputs.
- Full spatial statistics, atlas selection, or external resource provenance beyond the cited references.

## Upgrade path

- Swap in a pinned public benchmark slice or project dataset while preserving the deliverable names and schemas.
- Replace centroid-level fitting with a real Tangram mapping run and keep these QC files as lightweight summaries.
- Expand the spatial coherence block toward Squidpy graph metrics once the environment can support the heavier dependencies.
