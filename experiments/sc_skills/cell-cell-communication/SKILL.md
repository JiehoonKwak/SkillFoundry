---
name: cell-cell-communication-portable-skill
description: Use this experiment skill to run an end-to-end communication workflow that combines ligand-receptor scoring, group summaries, and optional spatial validation with public tools and explicit file contracts that work for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Run an end-to-end communication workflow that combines ligand-receptor scoring, group summaries, and optional spatial validation.

## Source adaptation

This experiment skill is adapted from `source_reference.md`, but it strips private Spatial Agent dependencies and keeps only shell-friendly public workflow steps.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every starter run reproducible from local files with no network access required at execution time.

## Focus terms

`cell-cell communication, ligand receptor, spatial`

## Default methods

- LIANA aggregation
- CellPhoneDB statistics
- spatial cross-check

## Recommended resource stack

- `LIANA API documentation`
- `LIANA+ paper`
- `CellPhoneDB methods documentation`
- `CellPhoneDB protocol paper`
- `Squidpy paper`
- `Squidpy neighborhood example`

## Core workflow

1. Validate that the toy or project input includes cell-level group labels, a ligand-receptor catalog, and an optional adjacency graph.
2. Compute group-level means and expression fractions before any interaction ranking.
3. Score each directed group pair with a lightweight LIANA-like aggregate from `min()` and geometric-mean components.
4. Run a tiny deterministic CellPhoneDB-like shuffle null to rank cell-type specificity.
5. Add a spatial adjacency cross-check when an edge list is available.
6. Keep the ranked interaction table, the smaller priority table, and the narrative report as separate deliverables.

## Deliverables

- `communication_results.tsv`
- `priority_pairs.tsv`
- `interpretation_report.md`

## QC artifacts

- `qc_group_expression.tsv`
- `qc_pair_nulls.tsv`
- `qc_spatial_support.tsv`
- `run_summary.json`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract and QC artifacts.
- `examples/toy_input.json` stores only raw synthetic cells, the ligand-receptor catalog, adjacency, and expected invariants.
- `scripts/run_cell_cell_communication.py` computes the starter outputs from scratch with deterministic local operations.
- `scripts/validate_outputs.py` verifies the deliverables and QC artifacts against `metadata.yaml`.
- `tests/test_cell_cell_communication.py` asserts numeric and structural invariants of the starter computation.
- `assets/README.md` records the starter boundary and upgrade path.

## Toy validation

- Run `python3 scripts/run_cell_cell_communication.py --input examples/toy_input.json --outdir scratch/cell-cell-communication_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/cell-cell-communication_toy_run` after the starter run.
- Run `python3 -m unittest discover -s tests -p 'test_*.py'` for the skill-local regression checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable names.
- Pin the database snapshot, orthology assumptions, preprocessing choices, and significance thresholds in the report or provenance notes.
- Escalate to real LIANA, CellPhoneDB, and Squidpy workflows when the dataset is large enough or when permutation depth and curated resources matter biologically.

## Starter scope

Computed: the starter truly computes group means, expression fractions, LIANA-like rank components, a tiny CellPhoneDB-like label-shuffle null, and adjacency support on the provided toy graph.
Approximated: LIANA is reduced to a two-component rank aggregate and CellPhoneDB is reduced to a small deterministic empirical null without database downloads or FDR correction.
Surrogate: the ligand-receptor catalog, cell matrix, and spatial graph remain synthetic starter inputs and are not a substitute for a real atlas-backed biological analysis.

## Failure handling

- Stop early when cell IDs are duplicated, ligand or receptor genes are missing from the input profiles, or adjacency references unknown cells.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve uncertainty in the report instead of inventing support that the starter did not compute.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance in `run_summary.json` for another shell-capable agent to rerun the starter exactly.
