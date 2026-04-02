---
name: cellphonedb-analysis-portable-skill
description: Use this experiment skill to run a portable CellPhoneDB-style ligand-receptor starter that computes statistical support, DEG relevance, and ranked interaction summaries with explicit file contracts for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Run CellPhoneDB in a portable way and summarize statistically supported interaction programs.

## Source adaptation

This experiment skill is adapted from `source_reference.md`, but it strips private Spatial Agent dependencies and keeps only shell-friendly public workflow steps.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every starter run reproducible from local files with no network access required at execution time.

## Focus terms

`cellphonedb, ligand receptor, single-cell`

## Default methods

- CellPhoneDB statistical mode
- DEG mode
- interaction ranking

## Recommended resource stack

- `CellPhoneDB methods and results documentation`
- `CellPhoneDB GitHub repository`
- `CellPhoneDB protocol paper`
- `cellphonedb-data repository`

## Core workflow

1. Validate that the input includes cell-level group labels, a ligand-receptor catalog, and a small DEG table.
2. Compute group-level means and expression fractions before any interaction testing.
3. Score each directed group pair with CellPhoneDB-like mean expression plus `min()` and geometric-mean summaries.
4. Run a tiny deterministic label-shuffle null to approximate CellPhoneDB statistical mode.
5. Mark DEG-relevant interactions and optional adjacency support as separate machine-readable QC layers.
6. Keep the significant means table, ranked interaction table, and narrative report as distinct deliverables.

## Deliverables

- `cellphonedb_significant_means.tsv`
- `interaction_ranked.tsv`
- `cellphonedb_report.md`

## QC artifacts

- `qc_group_expression.tsv`
- `qc_statistical_null.tsv`
- `qc_deg_support.tsv`
- `qc_adjacency_support.tsv`
- `run_summary.json`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract and QC artifacts.
- `examples/toy_input.json` stores only raw synthetic cells, the ligand-receptor catalog, DEG hints, adjacency, and expected invariants.
- `scripts/run_cellphonedb_analysis.py` computes the starter outputs from scratch with deterministic local operations.
- `scripts/validate_outputs.py` verifies the deliverables and QC artifacts against `metadata.yaml`.
- `tests/test_cellphonedb_analysis.py` asserts numeric and structural invariants of the starter computation.
- `assets/README.md` records the starter boundary and upgrade path.

## Toy validation

- Run `python3 scripts/run_cellphonedb_analysis.py --input examples/toy_input.json --outdir scratch/cellphonedb-analysis_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/cellphonedb-analysis_toy_run` after the starter run.
- Run `python3 -m unittest discover -s tests -p 'test_*.py'` for the skill-local regression checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable names and schemas.
- Pin the CellPhoneDB release, database snapshot, orthology assumptions, preprocessing choices, and significance thresholds in the report or provenance notes.
- Escalate to the real CellPhoneDB package when curated complexes, production permutation depth, or full scoring logic matter biologically.

## Starter scope

Computed: the starter truly computes group means, expression fractions, CellPhoneDB-like mean-expression tests, a tiny deterministic shuffle null, DEG relevance flags, and adjacency support on the provided toy graph.
Approximated: the statistical mode uses a bounded toy permutation count, DEG mode is reduced to a simple relevance table, and interaction ranking uses a lightweight heuristic instead of the full CellPhoneDB v5 scoring implementation.
Surrogate: the ligand-receptor catalog, DEG calls, and neighborhood graph remain synthetic starter inputs and are not a substitute for a real database-backed CellPhoneDB analysis.

## Failure handling

- Stop early when cell IDs are duplicated, ligand or receptor genes are missing from the input profiles, DEG annotations reference unknown groups, or adjacency references unknown cells.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve uncertainty in the report instead of inventing significance that the starter did not compute.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance in `run_summary.json` for another shell-capable agent to rerun the starter exactly.
