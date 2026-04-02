---
name: panel-design-portable-skill
description: Use this experiment skill to design a balanced marker panel for spatial transcriptomics or targeted assays using explicit evidence aggregation, redundancy control, and coverage balancing with public-tool-friendly file contracts for Codex, Claude Code, and similar shell-capable agents.
---
## Purpose

Design a balanced marker panel for spatial transcriptomics or targeted assays using explicit evidence and redundancy control.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every toy or pseudo-data run reproducible from local files.

## Focus terms

`panel design, marker genes, spatial, targeted assay`

## Default methods

- marker aggregation
- redundancy pruning
- coverage balancing

## Recommended resource stack

- `scGIST` paper for panel-size-constrained selection with prioritized genes
- `gpsFISH` paper for platform-aware targeted spatial panel selection
- `NS-Forest` docs or paper for minimum marker combinations and binary specificity logic
- `geneBasis` paper for manifold-preserving gene subset selection
- `PanglaoDB` and `CellTypist` for public marker priors
- `10x Xenium` panel-design docs for platform-facing follow-up constraints

## Core workflow

1. Aggregate marker evidence from literature, atlas, and database priors into one comparable score table.
2. Build a redundancy graph from near-duplicate marker profiles and keep one representative per component.
3. Satisfy per-cell-type quotas first, then spend the remaining slots on genes that improve undercovered compartments.
4. Leave behind machine-readable QC for aggregation, pruning, and balancing rather than only a prose report.

## Deliverables

- `panel_candidates.tsv`
- `panel_rationale.md`
- `platform_notes.md`

## Package scaffold

- `metadata.yaml` defines the deliverable contract for `panel-design`.
- `examples/toy_input.json` contains raw synthetic marker evidence, cell-type signal profiles, quotas, and expected invariants.
- `scripts/run_panel_design.py` computes weighted aggregation, redundancy pruning, quota selection, and coverage-balancing QC.
- `scripts/validate_outputs.py` checks the declared deliverables and QC artifacts.
- `tests/test_panel_design.py` asserts numeric invariants of the local starter path.
- `assets/README.md` explains what the starter really computes and where it remains a surrogate.

## Toy validation

- Run `python3 scripts/run_panel_design.py --input examples/toy_input.json --outdir scratch/panel-design_toy_run`.
- Run `python3 scripts/run_exercise.py --outdir scratch/panel-design_toy_run` to exercise the runner plus validator together.
- Run `python3 -m unittest discover -s tests -p 'test_*.py'` for the local contract and invariant checks.

## Real-run expectations

- Replace the toy evidence table with tissue-matched marker scores, atlas-derived prevalence, and platform-specific feasibility filters.
- Keep the same deliverable names so downstream validation still works.
- Record the exact atlas version, marker databases, platform chemistry, and panel size budget in the reports.

## Starter scope

- Truly computed: weighted marker aggregation from raw synthetic evidence, redundancy components from marker-profile similarity, and quota-plus-balance panel selection with QC traces.
- Approximated: probe feasibility and crowding are represented by lightweight scalar heuristics instead of sequence-level design or platform simulation.
- Surrogate only: no real atlas query, no classifier retraining, no direct probe/primer design, and no vendor panel optimizer.

## Minimum validation gates

- All declared deliverables exist and remain machine-readable.
- `panel_candidates.tsv` contains the required columns from `metadata.yaml`.
- `aggregated_marker_qc.tsv`, `redundancy_qc.tsv`, and `coverage_balance_qc.tsv` are written for every run.
- The tests confirm selected genes, redundancy outcomes, and final coverage totals on the toy input.

## Failure handling

- Stop early if quotas exceed the requested panel size or if no retained candidate can satisfy a required compartment.
- Preserve ambiguity in broad markers instead of pretending they are perfectly cell-type-specific.
- Treat the toy outputs as starter-path validation, not biological evidence for a real experiment.

## Hand-off contract

Before finishing, leave behind the declared deliverables, QC artifacts, and `run_summary.json` so another shell-capable agent can rerun or extend the panel-design pass.
