---
name: ligand-receptor-discovery-portable-skill
description: Use this experiment skill to discover biologically plausible ligand-receptor pairs and prioritize them with lightweight expression, interaction, and spatial evidence using a shell-friendly starter contract for Codex, Claude Code, and similar agents.
---
## Purpose

Discover biologically plausible ligand-receptor pairs and prioritize them with expression, interaction, and spatial evidence.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes private Spatial Agent tooling and keeps only portable shell-friendly workflow steps.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Prefer public package CLIs, Python APIs, and explicit file contracts.
- Keep every starter run reproducible from local files with no network access required at execution time.

## Focus terms

`ligand receptor discovery, communication, spatial`

## Default methods

- knowledge-base lookup
- expression filtering
- spatial validation

## Recommended resource stack

- `LIANA documentation`
- `LIANA+ paper`
- `CellPhoneDB v5 protocol`
- `CellPhoneDB analysis methods documentation`
- `OmniPath intercell network documentation`
- `OmniPath paper`
- `Squidpy neighborhood enrichment tutorial`
- `Squidpy paper`

## Core workflow

1. Validate that the input contains a raw cell matrix, group labels, a ligand-receptor catalog, and optional adjacency edges.
2. Compute group-level means and expression fractions before any ranking.
3. Score directed ligand-receptor candidates with expression filters plus lightweight knowledge-base support from the toy catalog.
4. Run a tiny deterministic label-shuffle null to estimate whether the observed expression score outranks simple relabelings.
5. Add adjacency support as a separate spatial QC layer and write final deliverables plus machine-readable intermediate QC.

## Deliverables

- `candidate_pairs.tsv`
- `evidence_table.tsv`
- `discovery_summary.md`

## QC artifacts

- `qc_group_expression.tsv`
- `qc_pair_scores.tsv`
- `qc_shuffle_null.tsv`
- `qc_spatial_support.tsv`
- `run_summary.json`

## Package scaffold

- `metadata.yaml` defines the stable deliverable contract and QC artifacts.
- `examples/toy_input.json` stores only raw synthetic cells, the ligand-receptor catalog, adjacency, resource aliases, and expected invariants.
- `scripts/run_ligand_receptor_discovery.py` computes the starter outputs from scratch with deterministic local operations.
- `scripts/validate_outputs.py` verifies deliverables and QC artifacts against `metadata.yaml`.
- `tests/test_ligand_receptor_discovery.py` asserts numeric and structural invariants of the starter computation.
- `assets/README.md` records the starter boundary and upgrade path.

## Toy validation

- Run `python3 scripts/run_ligand_receptor_discovery.py --input examples/toy_input.json --outdir scratch/ligand-receptor-discovery_toy_run` from this skill directory or with an absolute `--outdir`.
- Run `python3 scripts/validate_outputs.py --outdir scratch/ligand-receptor-discovery_toy_run` after the starter run.
- Run `python3 -m unittest discover -s tests -p 'test_*.py'` for the skill-local regression checks.

## Real-run expectations

- Replace `examples/toy_input.json` with project or public inputs that preserve the same deliverable names and schemas.
- Pin the exact resource snapshot, preprocessing choices, expression thresholds, and spatial graph construction in the report or provenance notes.
- Escalate to LIANA, CellPhoneDB, OmniPath-backed catalogs, and Squidpy or LIANA+ spatial workflows when richer databases, heteromeric complexes, or formal spatial statistics matter biologically.

## Starter scope

Computed: the starter truly computes group means, expression fractions, lightweight knowledge-base support from the toy catalog, deterministic ligand-receptor ranking, a tiny shuffle null, and adjacency support on the provided graph.
Approximated: the scoring is a small LIANA or CellPhoneDB-like surrogate built from min and geometric-mean expression plus resource weights, the null model is bounded to a tiny deterministic shuffle count, and the spatial check is an adjacency fraction rather than a full neighborhood-enrichment statistic.
Surrogate: the ligand-receptor catalog, expression matrix, and graph remain synthetic starter inputs and are not a substitute for a real biological discovery workflow.

## Failure handling

- Stop early when cell IDs are duplicated, ligand or receptor genes are missing from the toy matrix, or adjacency references unknown cells.
- Treat toy outputs as contract validation, not biological evidence.
- Preserve uncertainty in the report instead of inventing significance that the starter did not compute.

## Hand-off contract

Before finishing, leave behind the declared deliverables, the QC artifacts, and enough provenance in `run_summary.json` for another shell-capable agent to rerun the starter exactly.
