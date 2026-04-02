---
name: liana-analysis-portable-skill
description: Use this experiment skill to run a portable LIANA-style communication starter with explicit method selection, metadata harmonization, resource overlap summaries, and interpretable ranked interaction outputs for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Run LIANA or LIANA+ in a portable starter form and summarize prioritized ligand-receptor programs with explicit file contracts.

## Source adaptation

This experiment skill is adapted from `source_reference.md`, but it strips private Spatial Agent dependencies and keeps only shell-friendly public workflow steps.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every starter run reproducible from local files with no network access required at execution time.

## Focus terms

`liana, cell communication, ligand receptor`

## Default methods

- LIANA ranking
- resource harmonization
- condition-aware summaries

## Recommended resource stack

- `liana-py documentation`
- `LIANA GitHub repository`
- `LIANA+ paper`
- `LIANA comparison paper`
- `OmniPath paper`
- `CellPhoneDB statistical results documentation`

## Core workflow

1. Harmonize raw group labels, condition labels, and resource names before any scoring.
2. Compute group-level means and expression fractions from the raw toy matrix.
3. Score each directed group pair with explicit selected methods: `min_expr`, `geometric_mean`, `specificity`, and `shuffle_pvalue`.
4. Aggregate the selected method ranks with a deterministic mean-rank starter, then keep condition deltas and adjacency support as separate QC layers.
5. Write the ranked interaction table, resource overlap table, report, and machine-readable QC artifacts.

## Deliverables

- `liana_rankings.tsv`
- `resource_overlap.tsv`
- `liana_report.md`

## QC artifacts

- `qc_metadata_harmonization.tsv`
- `qc_group_expression.tsv`
- `qc_method_scores.tsv`
- `qc_condition_summary.tsv`
- `qc_adjacency_support.tsv`
- `run_summary.json`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract and QC artifacts.
- `examples/toy_input.json` stores only raw synthetic cells, the ligand-receptor catalog, harmonization aliases, adjacency, and expected invariants.
- `scripts/run_liana_analysis.py` computes the starter outputs from scratch with deterministic local operations.
- `scripts/validate_outputs.py` verifies deliverables and QC artifacts against `metadata.yaml`.
- `tests/test_liana_analysis.py` asserts numeric and structural invariants of the starter computation.
- `assets/README.md` records the starter boundary and upgrade path.

## Toy validation

- Run `python3 scripts/run_liana_analysis.py --input examples/toy_input.json --outdir scratch/liana-analysis_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/liana-analysis_toy_run` after the starter run.
- Run `python3 -m unittest discover -s tests -p 'test_*.py'` for the skill-local regression checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable names and schemas.
- Pin the exact LIANA release, resource snapshot, preprocessing choices, and selected methods in the report or provenance notes.
- Escalate to the real LIANA or LIANA+ stack when consensus over many back-end methods, spatial statistics, or multi-sample tensor workflows matter biologically.

## Starter scope

Computed: the starter truly computes metadata harmonization, group means, expression fractions, LIANA-like mean-rank prioritization across selected methods, a tiny deterministic shuffle null, resource overlap summaries, condition deltas, and adjacency support on the provided toy graph.
Approximated: the aggregation is a lightweight mean-rank surrogate rather than the full LIANA consensus implementation, the null model is bounded to a tiny permutation count, and specificity is a simple enrichment heuristic rather than every public scoring back-end.
Surrogate: the ligand-receptor catalog, canonical resources, and graph remain synthetic starter inputs and are not a substitute for a real OmniPath-backed or atlas-backed communication analysis.

## Failure handling

- Stop early when cell IDs are duplicated, ligand or receptor genes are missing from the toy matrix, harmonization aliases are incomplete, or adjacency references unknown cells.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve uncertainty in the report instead of inventing significance that the starter did not compute.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance in `run_summary.json` for another shell-capable agent to rerun the starter exactly.
