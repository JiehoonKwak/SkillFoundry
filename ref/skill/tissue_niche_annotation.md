# Tissue Niche Annotation for Single-Cell-Resolution Spatial Data

Annotate tissue niches or microenvironments in single-cell-resolution spatial data by combining cell-type labels, spatial neighborhoods, and local composition patterns.

## When to Use

- You already have per-cell observations with x/y coordinates.
- Cell types or broad cell states are already available, or can be transferred from a reference atlas.
- The goal is to assign each cell to a higher-level tissue niche such as immune aggregate, epithelial core, stromal boundary, or tumor-immune interface.
- You need a portable workflow that a general coding agent can run without relying on private Spatial Agent tools.

## Core Idea

Tissue niches are not the same as cell types. A niche label should reflect the **local cellular context** around a cell, not just the cell's own expression profile.

A practical portable workflow should:

1. Build a spatial neighbor graph from cell coordinates.
2. Summarize each cell's neighborhood composition by nearby cell types or states.
3. Optionally add compact expression or marker context for ambiguity resolution.
4. Group similar neighborhood profiles into candidate niches.
5. Interpret or annotate those niches using dominant cell types, boundary behavior, marker support, and spatial coherence.

## Public Method Anchors

- **Squidpy**
  - Spatial neighbor graph
  - Neighborhood enrichment
  - Co-occurrence and graph-aware spatial statistics

- **CellCharter**
  - Spatial clustering for cellular niches
  - Cluster characterization and neighborhood enrichment

- **Giotto / best-practices neighborhood analysis**
  - Local cell-type composition and spatial interaction analysis

## Portable Workflow

1. **Input checks**
   - Confirm single-cell resolution observations.
   - Verify `cell_id`, coordinates, and at least one label column such as `cell_type`.

2. **Build spatial graph**
   - Construct a deterministic kNN or radius graph from x/y coordinates.
   - Save per-cell neighbors and edge weights.

3. **Compute neighborhood composition**
   - For each cell, compute fractions or weighted counts of nearby cell types.
   - Derive boundary metrics such as heterogeneity, interface score, or dominant-neighbor margin.

4. **Assign candidate niche labels**
   - Use graph-smoothed composition signatures, compact clustering, or prototype scoring.
   - Keep the assigned niche plus a confidence or margin.

5. **Annotate the niche meaning**
   - Summarize dominant cell types and marker context for each niche.
   - Mark niche states like `immune_aggregate`, `stromal_boundary`, `epithelial_core`, or `mixed_interface` only when the evidence supports them.

6. **Validate**
   - Check that niche labels are spatially coherent.
   - Check that niche annotations reflect neighborhood composition rather than merely echoing the source cell type.
   - Flag ambiguous or low-support assignments explicitly.

## Expected Deliverables

- `niche_labels.tsv`
  - Per-cell niche assignments, confidence, dominant neighbor class, and boundary signals.

- `niche_markers.tsv`
  - Per-niche summary of enriched cell types, marker context, or composition statistics.

- `tissue_niche_report.md`
  - Narrative summary of graph construction, niche definitions, ambiguities, and caveats.

## Key Validation Questions

- Does the workflow derive niche labels from neighborhood context instead of simply copying `cell_type`?
- Are boundary or interface cells identifiable from mixed neighborhoods?
- Are niche assignments reproducible on deterministic toy inputs?
- Are ambiguity and low-support cases preserved instead of overclaimed?
