# Cell-Cell Communication Decomposition Notes

This file keeps only a sanitized decomposition of the original task.
It does not depend on private Spatial Agent tools, hidden wrappers, or notebook-only automation.

## Public workflow decomposition

1. Inspect cell-level inputs and confirm that every cell has a group label, an expression profile, and optional spatial coordinates.
2. Summarize group-level expression for each ligand and receptor candidate.
3. Score each source-group to target-group ligand-receptor pair with simple deterministic statistics:
   - a LIANA-like magnitude view based on `min(ligand_expr, receptor_expr)`,
   - a LIANA-like magnitude view based on geometric mean,
   - a CellPhoneDB-like null model based on shuffled group labels.
4. If an adjacency graph is available, compute an optional spatial cross-check over neighboring cell pairs.
5. Emit three stable deliverables for downstream agents:
   - `communication_results.tsv`
   - `priority_pairs.tsv`
   - `interpretation_report.md`
6. Emit machine-readable QC so another shell-capable agent can audit the starter run without re-reading the report text.

## Public replacements for the removed private stack

- `liana-py` documentation and paper for ligand-receptor ranking and rank aggregation.
- `CellPhoneDB` documentation and protocol paper for thresholded means plus permutation-based cell-type specificity.
- `Squidpy` documentation and paper for neighborhood-aware cross-checks and spatial context.

## Starter boundary

The experiment package intentionally uses tiny synthetic inputs and deterministic local computation.
It approximates the public methods closely enough to exercise the workflow contract, but it is not a substitute for a full LIANA, CellPhoneDB, or Squidpy analysis on real data.
