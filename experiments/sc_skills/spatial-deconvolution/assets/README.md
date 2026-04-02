# Assets

This directory belongs to the experiment skill `spatial-deconvolution`.

The current refine pass upgrades the package from a static contract to a tiny runnable starter. It still uses deterministic synthetic inputs only, so no biological claims should be drawn from these assets.

## What the starter computes

- Shared-gene alignment between the spot matrix and reference signatures.
- A small `cell2location`-like NNLS abundance estimate on library-normalized counts.
- A small `DestVI`-like surrogate path that smooths counts across a toy spatial neighbor graph before log-NNLS fitting.
- A uniform baseline comparison and machine-readable QC summaries.

## What the starter approximates

- `cell2location` is approximated with deterministic NNLS instead of Bayesian posterior inference with tissue priors.
- `DestVI` is approximated with log-scale stabilization plus coordinate-based smoothing instead of latent variable training and cell-state modeling.
- Any truth-aware metric comes from synthetic mixture tables supplied in the toy input, not from an external benchmark or hidden evaluation server.

## What remains outside scope

- Real atlas construction or reference signature learning.
- Posterior intervals, uncertainty calibration, and model selection on real tissue.
- GPU-backed training, hyperparameter tuning, and biological interpretation of spatial niches.

## Upgrade path

- Replace the toy input with a pinned public dataset slice while preserving the same deliverable names and schemas.
- Swap the surrogate NNLS paths for the real public methods when the runtime budget and environment permit it.
- Keep the QC files and summary sections explicit so agents can tell which parts are exact computations and which remain approximations.
