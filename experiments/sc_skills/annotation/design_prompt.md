Design or refine an experiment-only portable skill for the task: Cell Type Annotation for Single-Cell-Resolution Spatial Data.

Task goal:
Replace the current toy-centric annotation skill with an operator-grade workflow that a general coding agent can actually follow on real data. The package must still keep deterministic local validation, but the skill itself should now read like a concrete execution playbook rather than a toy-method summary.

Primary scenario:
- Single-cell-resolution spatial platforms such as MERFISH, Xenium, CosMx, and seqFISH.
- Not for spot-based platforms such as Visium, Slide-seq, ST, or 10x Spatial Gene Expression. Those must be routed to `spatial-deconvolution`.

High-level outcome:
- A shell-capable agent should be able to open `SKILL.md` and execute the workflow step by step:
  - detect platform type
  - explore the dataset
  - search for a CELLxGENE/Census reference
  - download or pin the reference
  - preprocess query + reference
  - run Harmony-style integration
  - transfer labels hierarchically
  - handle missing tools with explicit fallback/bootstrap instructions

Critical environment fact:
- The local default interpreter is Python 3.13.
- `cellxgene-census` is currently blocked on Python 3.13 because published wheels/metadata require Python `<3.13`.
- The redesigned skill must explain this concretely and provide exact environment-creation commands for a supported Python 3.10-3.12 runtime, plus an alternate metadata-export / manual-browser route when such an environment cannot be created locally.

Constraints:
- Only create or modify files under `experiments/sc_skills/annotation`.
- Do not edit `registry/`, `skills/`, `site/`, `README.md`, `experiments.md`, or planning files.
- Keep the result compatible with Codex, Claude Code, and similar shell-capable agents.
- Use `source_reference.md` only as decomposition guidance; remove private Spatial Agent tool assumptions.
- Preserve the declared deliverable names and schemas so repository-level validation still passes.
- Keep local validation deterministic and offline, but do not let toy references or toy catalogs dominate the real-run instructions.

Hard requirements for this redesign:
1. `SKILL.md` must become a concrete runbook, not just a high-level summary.
2. The runbook must include explicit step-by-step instructions with:
   - what to inspect
   - what commands or Python snippets to run
   - what files should be produced at each step
   - what decisions to make based on the observed results
3. The runbook must include explicit search guidance for reference discovery:
   - how to form queries
   - how to inspect CELLxGENE/Census dataset metadata
   - how to decide which reference to keep
4. The runbook must include explicit missing-tool handling:
   - how to check whether required packages are present
   - exact install/bootstrap suggestions
   - what bundled helper scripts to use when the heavy public tool is not available
   - what to do specifically when the host interpreter is Python 3.13 and `cellxgene-census` cannot be installed
5. The package must include at least one helper script for real-run preparation, such as:
   - runtime/tooling checks
   - reference-ranking from exported metadata
   - run-plan scaffolding
   Use whatever filenames fit best, but make them part of the package contract and document them.
6. The toy validation path should remain deterministic, but it should be clearly separated from the real-run workflow in the docs.
7. Tests must assert workflow behavior, not only file existence. Include at least one test for the spot-based rejection/reroute path and at least one test that the search/reference-selection logic is actually used.

Workflow content requirements:

1. Platform detection and dataset exploration
   - inspect shape, obs/var columns, coordinates, value range, candidate sample/batch columns
   - determine single-cell spatial vs spot-based
   - explain what metadata to capture before moving on

2. Query preprocessing
   - describe when to normalize and when to skip
   - show representative Scanpy commands or equivalent Python
   - record expected intermediate artifacts or QC tables

3. Reference discovery from CELLxGENE/Census
   - show concrete search logic such as `{species} {tissue} normal`
   - explain how to inspect dataset tables / metadata
   - explain what to prefer: healthy/normal, matching species, same tissue, good annotations, sufficient cells
   - explain when to reject a candidate reference

4. Reference download and preparation
   - include example code for `cellxgene_census.download_source_h5ad(...)`
   - pin dataset ID, version, and provenance
   - explain how to harmonize cell-type columns before transfer

5. Integration and label transfer
   - give concrete steps around gene intersection, normalization, PCA, Harmony
   - explain broad-label first, subtype second
   - explain how to review low-confidence calls

6. Missing-tool fallback
   - if `cellxgene_census` is missing, explain how to bootstrap it or fall back to a ranked metadata export path
   - if Harmony tooling is missing, explain how to use the bundled surrogate path and what accuracy limitations it has
   - include a supported-environment recipe using `conda` or another explicit Python 3.10-3.12 environment creator
   - include a browser/manual CELLxGENE Discover route if the Python package path is blocked

7. Deliverables and stop conditions
   - explicitly list required outputs
   - explain when to stop because the reference is poor, genes do not overlap, coordinates are missing, or labels are not trustworthy
   - explicitly route tissue niche analysis to `tissue-niche-annotation` instead of doing it here

Preferred public anchors:
- CELLxGENE Census quick start
- CELLxGENE datasets table tutorial
- `cellxgene_census.download_source_h5ad`
- Scanpy preprocessing
- `scanpy.external.pp.harmony_integrate`
- CellTypist and/or Azimuth as label-transfer design anchors
- `cellxgene-census` PyPI package metadata / install constraints
- CELLxGENE Discover web dataset browser as the manual-search fallback

Package-level expectations:
- Keep or improve the existing runnable toy starter, but make the docs much more actionable for real data.
- Update `refs.md` with any additional primary sources used.
- Update `assets/README.md` to explain the real-run vs surrogate boundary clearly.
- Refresh `metadata.yaml`, `examples/`, `scripts/`, and `tests/` so they match the new operator-grade workflow.
- Keep this skill strictly scoped to cell-type annotation. Do not perform tissue niche annotation in this package; instead document the handoff to `experiments/sc_skills/tissue-niche-annotation`.
