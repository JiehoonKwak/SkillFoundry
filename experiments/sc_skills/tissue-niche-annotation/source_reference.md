# Tissue Niche Annotation for Single-Cell-Resolution Spatial Data

Annotate tissue niches or microenvironments in single-cell-resolution spatial data by combining cell-type labels, spatial neighborhoods, and local composition patterns.

## Use case boundary

- Use this workflow only when each observation is a single cell with spatial coordinates.
- The workflow assumes an existing cell-type or state annotation column.
- If labels are missing or low-confidence, route back to the annotation workflow first.
- If the platform is spot-based, use spatial deconvolution or spot-domain workflows instead.

## Public method anchors

- Squidpy for spatial neighbor graphs and neighborhood enrichment.
- CellCharter for neighborhood aggregation and clustering.
- UTAG for graph-aware microenvironment grouping.
- Single-cell spatial best-practices material for neighborhood analysis and interpretation.

## Operator checklist

1. Inspect the AnnData object and capture:
   - coordinate keys
   - label column
   - sample or slide key
   - platform hints
   - whether those labels came from `annotation_table.tsv` and whether `marker_status` suggests they are trustworthy enough for downstream niche work
2. Run the runtime checker and record missing dependencies.
3. Generate a run plan before graph construction so graph parameters and expected outputs are explicit.
4. Build a spatial graph with Squidpy when available; otherwise use the bundled surrogate path only as provisional fallback.
5. Compute neighborhood composition and, when possible, neighborhood enrichment.
6. Derive candidate niches using CellCharter, UTAG, or the documented fallback path.
7. Summarize niche meaning with dominant cell types, boundary behavior, marker context, and ambiguity notes.

## Required intermediate artifacts

- `00_runtime_check.json`
- `01_dataset_profile.json`
- `02_niche_run_plan.json`
- `03_spatial_graph.tsv` or a graph-ready H5AD
- `04_neighborhood_profiles.tsv`
- `05_niche_candidates.tsv`
- `06_niche_review.md`

## Key decisions

- Prefer `kNN` when cell density is uneven or coordinate units are uncertain.
- Prefer `radius` only when the coordinate system is calibrated and the biological scale is known.
- Keep interface or boundary cells explicit; do not force every cell into a confident core niche.
- A niche label must reflect local context rather than simply repeating the focal cell's label.

## Missing-tool behavior

- If `squidpy` is missing, bootstrap it when possible; otherwise keep the graph stage provisional and use the bundled deterministic fallback.
- If `cellcharter` is missing, use neighborhood profiles plus heuristic or prototype-based grouping and mark the result accordingly.
- If `utag` is missing, use the fallback graph-aware grouping path and record that the run is not a formal UTAG execution.
- If the upstream annotation is weak, route back to the `annotation` skill instead of trying to repair label transfer here.

## Deliverables

- `niche_labels.tsv`
- `niche_markers.tsv`
- `tissue_niche_report.md`
