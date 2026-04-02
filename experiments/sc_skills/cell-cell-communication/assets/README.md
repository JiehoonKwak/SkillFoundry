# Assets

This directory belongs to the experiment skill `cell-cell-communication`.

The current refine pass upgrades the package from a contract-only stub to a tiny runnable starter with deterministic local computation on synthetic cells, a small ligand-receptor catalog, and a toy adjacency graph.

## What the starter computes

- Group-level mean expression and expression fractions for every gene present in the toy input.
- A lightweight LIANA-like aggregate that combines `min(ligand_mean, receptor_mean)` and geometric-mean ranks.
- A tiny CellPhoneDB-like label-shuffle null ranking for cell-type specificity.
- An adjacency-based spatial support fraction over the provided toy neighbor graph.

## What is approximated

- The LIANA consensus is reduced to two magnitude-oriented components instead of the full method panel and resource harmonization stack.
- The CellPhoneDB null model uses only a small deterministic set of shuffles and does not reproduce the full database, multi-subunit handling, or multiple-testing workflow.
- The spatial cross-check is a simple graph support calculation, not a full Squidpy neighborhood enrichment or ligand-receptor permutation analysis.

## What the starter does not validate

- Biological correctness on real single-cell or spatial transcriptomics data.
- Curated database completeness, orthology conversion, or complex-subunit logic.
- Differential-condition modeling, atlas integration, or large-scale permutation calibration.

## Upgrade path

- Replace the toy cells with a pinned real or benchmark dataset slice while keeping the same deliverable filenames.
- Swap the synthetic ligand-receptor catalog for a curated public resource such as OmniPath or CellPhoneDB.
- Promote the spatial cross-check to a real Squidpy neighborhood or ligand-receptor analysis once the runtime budget and environment allow it.
