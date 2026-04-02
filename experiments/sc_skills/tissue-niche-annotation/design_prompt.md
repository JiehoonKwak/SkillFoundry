Design or refine an experiment-only portable skill for the task: Tissue Niche Annotation for Single-Cell-Resolution Spatial Data.

Task goal:
Turn this package into an operator-grade workflow that a general coding agent can execute on real single-cell-resolution spatial data. The skill must be highly detailed, step-by-step actionable, and explicit about which public tools to use, how to inspect inputs, how to save intermediate artifacts, and how to recover if required tooling is missing.

Primary scenario:
- Single-cell-resolution spatial data with per-cell coordinates and an existing cell-type or state annotation column.
- Typical platforms: MERFISH, Xenium, CosMx, seqFISH.
- This skill assumes annotation is already available. If cell types are missing or low-confidence, the workflow must explicitly tell the agent to return to the cell-type annotation workflow first.
- This skill must remain separate from cell-type annotation. Do not redesign it into a reference-search or label-transfer workflow.

Required outcome:
- `SKILL.md` must become a concrete real-run playbook, not a short toy-method summary.
- The package must keep deterministic local validation, but the user-facing workflow must be real-workflow-first.
- The package must include helper surfaces for runtime checks and niche-run preparation so the workflow is executable rather than only descriptive.

Constraints:
- Only create or modify files under `experiments/sc_skills/tissue-niche-annotation`.
- Do not edit `registry/`, `skills/`, `site/`, `README.md`, `experiments.md`, or planning files.
- Keep the result compatible with Codex, Claude Code, and similar shell-capable agents.
- Preserve declared deliverable names and schemas so repository-level validation still passes.
- Keep toy validation deterministic, but do not let synthetic niche prototypes dominate the real-run instructions.

Hard requirements:
1. `SKILL.md` must contain a numbered real-run workflow. For each step, include:
   - purpose
   - exact commands or Python snippets
   - expected intermediate outputs
   - decision criteria and stop conditions
2. The workflow must explicitly show how to inspect the input AnnData and determine whether the data are ready for niche analysis.
3. The workflow must explicitly show how to choose and build the spatial graph:
   - kNN vs radius
   - per-sample handling
   - what QC to save
4. The workflow must explicitly show how to compute neighborhood features and derive candidate niches using public-tool anchors such as Squidpy, CellCharter, or UTAG.
5. The workflow must explicitly explain how to interpret niches using local composition, boundary behavior, and marker context, and how to avoid simply renaming cell types as niches.
6. The package must include at least two helper scripts beyond the toy runner, for example:
   - a runtime/tooling checker
   - a niche-run planner or graph-parameter helper
   Use whatever filenames fit best, but make them part of the package contract and document them.
7. Missing-tool fallback is mandatory:
   - if Squidpy is missing, explain how to install it and how to use bundled local graph utilities instead
   - if CellCharter or UTAG are missing, explain how to use the bundled deterministic fallback scoring path and what limitations it has
   - if upstream annotation is inadequate, explicitly route back to the annotation workflow
   - if a dedicated environment is required, provide exact commands to create it
8. Tests must validate workflow behavior, not only file existence. Include assertions for:
   - dependence on neighborhood context rather than copied cell types
   - presence of operator-grade real-run instructions in docs
   - helper-script behavior on deterministic input

Detailed workflow content requirements:

1. Dataset triage
   - inspect shape, coordinates, sample IDs, label columns, and possible batch columns
   - state the minimum metadata required before niche analysis can start
   - explicitly reject data without coordinates or without trustworthy cell-type/state labels

2. Spatial graph construction
   - show representative Squidpy or Python code
   - explain when to prefer kNN vs radius
   - explain how to handle multiple samples or slides
   - save graph/QC outputs

3. Neighborhood feature construction
   - compute local cell-type composition
   - optionally smooth over the graph
   - explain when to prefer raw vs smoothed composition
   - save intermediate neighborhood tables

4. Candidate niche discovery
   - use CellCharter-style aggregation, UTAG-style grouping, or a documented deterministic surrogate
   - explain how to choose niche number or prototype structure
   - keep interface or boundary cells explicit

5. Niche interpretation
   - summarize dominant cell types, marker context, interface behavior, and ambiguity
   - explain what evidence is needed before calling a niche `immune_aggregate`, `stromal_boundary`, `epithelial_core`, `mixed_interface`, etc.

6. Tool recovery / bootstrap
   - show concrete package checks
   - give concrete install commands or environment suggestions
   - if the public tool cannot be installed, point to the bundled helper scripts or surrogate path

7. Deliverables and stop conditions
   - explicitly list required outputs
   - explain when to stop because the annotation is weak, graph construction is unstable, or niche separation is not credible

Preferred public anchors:
- Squidpy `spatial_neighbors`
- Squidpy `nhood_enrichment`
- Squidpy tutorials for single-cell-resolution spatial data
- CellCharter docs and paper
- UTAG repository/paper
- single-cell spatial best-practices documentation

Package-level expectations:
- Keep or improve the deterministic toy starter and its tests.
- Make the docs and helper scripts clearly usable on real data.
- Update `refs.md` with any additional primary sources used.
- Update `assets/README.md` to explain the real-run vs surrogate boundary clearly.
- Refresh `metadata.yaml`, `examples/`, `scripts/`, and `tests/` so they match the stronger workflow.
- Make the dependency boundary between this package and `annotation` explicit: this package consumes cell-type/state labels and should not repeat reference search or label transfer.
