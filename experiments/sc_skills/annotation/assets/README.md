# Assets Boundary

This experiment package has two clearly separated operating modes.

## Real-run path

- `SKILL.md` is the operator-facing runbook for MERFISH, Xenium, CosMx, and seqFISH style single-cell spatial datasets.
- `scripts/check_runtime.py` records which public packages are present, whether the current Python can run `cellxgene-census`, and which bootstrap commands to use when it cannot.
- `scripts/plan_census_env.py` tells the agent whether to:
  - use the host Python directly,
  - create a dedicated `conda` Python 3.12 environment,
  - create a `python3.12` venv,
  - or switch to the Discover browser fallback.
- `scripts/rank_reference_metadata.py` ranks an exported CELLxGENE/Census metadata table offline and explains why each candidate is kept or rejected.
- The real-run instructions assume that the agent will inspect a real query dataset, search or export real Discover/Census metadata, pin a real reference dataset ID, and write step-specific QC artifacts before trusting transferred labels.

## Deterministic surrogate path

- `scripts/run_annotation.py` and `examples/toy_input.json` remain the offline validation path.
- The toy run still computes platform classification, query normalization, reference ranking, and hierarchical label transfer with marker review.
- The toy outputs are deterministic and are meant to validate workflow shape, output schemas, and fallback logic.

## What remains approximate

- The toy catalog is not a live CELLxGENE/Census catalog.
- The Harmony step is a deterministic batch-centering surrogate, not iterative Harmony integration through `scanpy.external.pp.harmony_integrate`.
- The label-transfer stage uses reference-neighbor voting plus marker gaps instead of a trained CellTypist or Azimuth reference model.

## What the toy path does not prove

- Biological correctness on a real tissue or disease context.
- Reference fitness for a real study, including ontology cleanup, donor balance, or disease-control comparability.
- Any tissue niche or microenvironment claim. Those belong in the separate `tissue-niche-annotation` skill.

## How to use the boundary

- Use the real-run path whenever you have a real `.h5ad` or equivalent single-cell spatial dataset and access to public references.
- On Python 3.13 hosts, prefer the supported Python 3.12 route recommended by `scripts/plan_census_env.py`; if that is impossible, use the Discover browser metadata-export fallback.
- Use the deterministic surrogate path to validate the package contract, rehearse the workflow, or keep moving when the environment cannot yet install the public reference-transfer stack.
- Do not present the surrogate outputs as final biological evidence. Treat them as a portability and decision-support scaffold only.
