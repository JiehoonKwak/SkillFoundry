# Assets

This directory belongs to the experiment skill `panel-design`.

The current starter uses synthetic marker evidence only. It is designed to exercise the panel-selection workflow locally and deterministically.

## What the starter computes

- Weighted aggregation of literature, atlas, and database support into a candidate marker score.
- Redundancy components from near-duplicate cell-type signal profiles.
- Quota-based anchor selection followed by a final coverage-balancing fill step.
- Machine-readable QC tables for aggregation, pruning, balancing, and a compact `run_summary.json`.

## What is approximated

- Probe suitability is represented by a toy `probe_design_score`.
- Panel crowding is represented by a scalar `crowding_penalty`.
- Evidence values are synthetic priors rather than live pulls from PanglaoDB, CellTypist, CELLxGENE, or vendor tools.

## What remains outside scope

- Sequence-level probe or primer design.
- Real atlas retrieval, cross-dataset harmonization, or classifier validation.
- Platform-specific chemistry simulation, imaging constraints, or vendor approval rules.

## Upgrade path

- Replace the toy evidence table with tissue-matched marker rankings and real platform feasibility filters.
- Keep the deliverable names stable so the portable contract still validates.
- Preserve the QC outputs when promoting the starter toward a stronger public-method wrapper.
