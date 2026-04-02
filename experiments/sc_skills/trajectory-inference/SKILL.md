---
name: trajectory-inference-portable-skill
description: Use this experiment skill to infer trajectories, pseudotime, RNA velocity, and fate probabilities with explicit preprocessing checkpoints and portable shell-first outputs for Codex, Claude Code, and similar agents.
---
## Purpose

Infer trajectories, pseudotime, RNA velocity, and fate probabilities with explicit preprocessing and interpretation checkpoints.

## Source adaptation

This experiment skill is derived from `source_reference.md`, but it removes Spatial Agent private tools, hidden wrappers, and notebook-only orchestration.

## Agent compatibility

- Compatible with Codex, Claude Code, and similar shell-capable agents.
- Uses deterministic local Python plus explicit TSV, JSON, and Markdown outputs.
- Avoids network access and heavy single-cell runtimes in the starter path.

## Default methods

- k-NN graph construction on toy coordinates
- root-to-cell shortest-path pseudotime
- neighbor-delta RNA velocity surrogate from synthetic spliced and unspliced programs
- absorbing Markov-chain fate probabilities toward two terminal states

## Deliverables

- `trajectory_coordinates.tsv`
- `fate_probabilities.tsv`
- `trajectory_report.md`

## Intermediate QC

The runner also emits:

- `qc_preprocessing.json`
- `qc_knn_graph.tsv`
- `qc_velocity_surrogate.tsv`
- `qc_transition_matrix.tsv`
- `run_summary.json`

## Starter scope

This starter truly computes a k-NN graph, shortest-path pseudotime, a directed transition matrix, and absorbing-state fate probabilities on the toy graph. It approximates RNA velocity with two synthetic spliced and unspliced programs plus neighbor-expression differences. It remains a surrogate for the full public scVelo and CellRank stack, so it does not perform moment estimation, dynamical model fitting, uncertainty propagation, or real AnnData preprocessing.

## Recommended resource stack

- `Scanpy DPT documentation`
- `PAGA paper`
- `scVelo documentation and paper`
- `CellRank tutorials and paper`

## Core workflow

1. Validate the root cell, terminal cells, and graph connectivity on the raw toy inputs.
2. Build the symmetric k-NN graph and compute shortest-path pseudotime from the root.
3. Derive a lightweight velocity direction from synthetic spliced and unspliced programs plus forward-neighbor expression deltas.
4. Convert forward transitions into a tiny absorbing Markov chain and summarize cell-level fate probabilities.
5. Write deliverables plus machine-readable QC artifacts for downstream checks.

## Example

```bash
python3 experiments/sc_skills/trajectory-inference/scripts/run_exercise.py \
  --outdir scratch/trajectory-inference_toy_run
```

## Validation

- Runner: `python3 experiments/sc_skills/trajectory-inference/scripts/run_trajectory_inference.py --input experiments/sc_skills/trajectory-inference/examples/toy_input.json --outdir scratch/trajectory-inference_toy_run`
- Validator: `python3 experiments/sc_skills/trajectory-inference/scripts/validate_outputs.py --outdir scratch/trajectory-inference_toy_run`
- Tests: `python3 -m unittest discover -s experiments/sc_skills/trajectory-inference/tests -p 'test_*.py'`

## Interpretation checkpoints

- Treat `qc_preprocessing.json` as the first gate: the graph must stay connected and the root and terminal cells must be present.
- Use `trajectory_coordinates.tsv` to confirm monotonic pseudotime along the intended trunk and branches before reading fate outputs.
- Use `qc_velocity_surrogate.tsv` only as a directional cue; it is not a substitute for full RNA-velocity parameter inference.
- Use `fate_probabilities.tsv` and `trajectory_report.md` for a compact branch summary and explicit caveats.

## Real-run expectations

- Replace the synthetic cells with project-specific coordinates or latent features while keeping the deliverable names unchanged.
- Swap the surrogate velocity step for a public scVelo or CellRank workflow when real spliced and unspliced counts are available.
- Record dataset provenance, preprocessing choices, and model parameters in the report before treating any result as biological evidence.
