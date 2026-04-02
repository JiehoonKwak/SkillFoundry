# Assets

This directory belongs to the experiment skill `trajectory-inference`.

## What the starter computes

- A symmetric k-NN graph on ten synthetic cells with one root, one branchpoint, and two terminal states.
- Root-to-cell shortest-path pseudotime on that graph.
- A lightweight RNA velocity surrogate from two synthetic spliced and unspliced programs plus forward-neighbor expression differences.
- Fate probabilities from a tiny absorbing Markov chain whose absorbing states are the two terminal cells.

## What is approximated

- The velocity step is a deterministic directional heuristic. It is not scVelo moment estimation, stochastic velocity, or a dynamical model fit.
- The fate step captures the CellRank-style idea of absorbing terminal states, but not full kernel composition, uncertainty propagation, or estimator fitting on real AnnData objects.

## What remains outside starter scope

- Real count-matrix preprocessing, highly variable gene selection, neighborhood construction in latent transcriptomic space, and uncertainty-aware velocity kernels.
- Any biological interpretation beyond the synthetic branch geometry and invariants encoded in `examples/toy_input.json`.

## Upgrade path

- Replace the synthetic inputs with a pinned public benchmark slice or project-specific latent coordinates and spliced/unspliced layers.
- Keep the same deliverable names while promoting the surrogate velocity and fate steps to real scVelo and CellRank runs.
- Preserve machine-readable QC outputs so agents can still verify preprocessing and interpretation checkpoints locally.
