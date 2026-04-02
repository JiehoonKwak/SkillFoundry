# Assets

This directory belongs to the experiment skill `ligand-receptor-discovery`.

The current refine pass upgrades the package from a contract-only toy into a tiny runnable starter over deterministic synthetic data. No real biological claims should be drawn from these assets.

## What the starter truly computes

- Group-level means and expression fractions from the raw toy cell matrix.
- Lightweight ligand-receptor support scores that combine expression summaries, a toy catalog resource lookup, a tiny label-shuffle null, and adjacency support.
- Machine-readable QC tables that expose the intermediate computations used to rank the final candidates.

## What remains approximated

- Resource support is derived from a tiny local catalog plus alias-and-weight rules instead of a real OmniPath or CellPhoneDB database snapshot.
- The null model uses only 31 deterministic shuffles and no multiple-testing correction.
- Spatial validation is reduced to adjacency support on a small graph rather than a full Squidpy neighborhood statistic or LIANA+ local metric.

## What the starter does not validate

- Biological correctness on real single-cell or spatial datasets.
- Heteromeric complexes, species translation, or atlas-specific annotation issues.
- Production-scale ranking stability, database versioning, or downstream pathway interpretation.

## Upgrade path

- Replace the toy input with a pinned public or project dataset that preserves the same deliverable names and schemas.
- Swap the local catalog surrogate for a real OmniPath, CellPhoneDB, or LIANA resource export and record the exact snapshot used.
- Replace adjacency-fraction support with a real spatial-neighborhood workflow when coordinates and runtime budget allow it.
