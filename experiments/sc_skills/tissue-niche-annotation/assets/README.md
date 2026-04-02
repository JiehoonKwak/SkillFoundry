# Assets

This directory belongs to the experiment skill `tissue-niche-annotation`.

## Real-run versus surrogate boundary

This package now presents a real-run-first workflow, but the bundled runnable path is still deterministic and synthetic so local validation stays stable.

## What the bundled starter truly computes

- a spatial kNN graph from raw toy coordinates
- per-cell neighborhood composition vectors
- one-step graph smoothing over those profiles
- prototype-based niche scores with explicit boundary metrics
- niche-vs-rest marker summaries on small synthetic inputs

## What remains approximate

- Squidpy neighborhood enrichment is represented by graph-aware local composition rather than permutation testing
- CellCharter-style aggregation and clustering are represented by explicit neighborhood features and niche prototype scoring
- UTAG-style grouping is represented by graph-aware deterministic surrogate logic rather than the full external package runtime

## What to do on real data

- Use `scripts/check_runtime.py` to verify the public runtime stack.
- Use `scripts/plan_niche_run.py` to scaffold graph parameters, expected outputs, and stop conditions.
- If `conda` is available, prefer a dedicated `python=3.12` environment for Squidpy, CellCharter, and UTAG rather than mixing them into an unrelated host environment.
- Start from a trustworthy label column, ideally `predicted_label` or `broad_label` produced by the separate `annotation` skill.
- Prefer public-tool execution with Squidpy, CellCharter, and UTAG when they are available.
- Use the bundled surrogate only when the required public tool cannot be installed, and label the result as provisional.
