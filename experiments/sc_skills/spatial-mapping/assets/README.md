# Assets

This directory belongs to the experiment skill `spatial-mapping`.

The current starter uses deterministic synthetic inputs only. No real biological claims should be drawn from these assets.

## Starter boundary

- Computed locally: gene intersection, dropped-feature QC, library-size normalization, `log1p`, reference-centroid cosine similarity, top-label margins, 3-NN niche majority voting, and marker consistency checks.
- Approximated: Tangram's full training loop, density priors, and projected-gene workflow are replaced by a centroid-level cosine surrogate so the starter stays portable and finishes in under 5 seconds.
- Surrogate: the marker heuristics, confidence thresholds, and toy reference profiles are synthetic, so the package validates workflow shape and QC logic rather than biological performance on a public benchmark.

## What the starter validates

- The runnable contract between `examples/`, `scripts/`, `metadata.yaml`, and the declared deliverables.
- Machine-readable QC artifacts that another shell-capable agent can inspect or test.
- Deterministic numerical behavior for gene intersection, mapping scores, review statuses, and neighborhood coherence on tiny toy data.

## What the starter does not validate

- Biological correctness on real spatial transcriptomics datasets.
- Tangram training dynamics, PyTorch optimization behavior, density priors, or projected-gene outputs.
- Full atlas selection, histology-aware neighborhood metrics, or external resource provenance beyond the cited references.

## Upgrade path

- Swap in a pinned public benchmark slice or project dataset while preserving the deliverable names and schemas.
- Replace centroid-level cosine scoring with a real Tangram mapping run and keep the QC files as lightweight review summaries.
- Expand the neighborhood review block toward richer Squidpy graph metrics once the environment can support the heavier dependencies.
