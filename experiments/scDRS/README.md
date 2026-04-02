# scDRS Experiment

This directory holds an experiment-only skill for running the scDRS workflow from GWAS summary statistics to single-cell disease-relevance scoring. It is intentionally isolated from the production skill registry and exists to prototype a portable, agent-friendly workflow under `experiments/scDRS/`.

## What is different from production skills

- It is not registered in `registry/` and is not part of `skills/`.
- It is written for generic coding agents such as Codex and Claude Code, so every step uses explicit files, CLI commands, checkpoints, and failure handling instead of private tool calls.
- It includes experiment-scoped references, a resource inventory, a deliverables checklist, and a fill-in execution template so the workflow can be validated before any future productionization decision.

## Directory contents

- `scdrs_gwas_to_singlecell_skill/SKILL.md`: main portable skill instructions
- `scdrs_gwas_to_singlecell_skill/refs.md`: cited sources and why they matter
- `scdrs_gwas_to_singlecell_skill/resource_inventory.md`: curated canonical materials
- `scdrs_gwas_to_singlecell_skill/deliverables_checklist.md`: expected runtime artifacts and validation gates
- `scdrs_gwas_to_singlecell_skill/templates/scdrs_execution_checklist.md`: reusable run template for a real execution

## Validation scope

The experiment is validated by a repository-level regression test that checks the required files and stable sections exist. It does not register the skill, alter site data, or claim runtime verification of MAGMA or scDRS in the current repository environment.
