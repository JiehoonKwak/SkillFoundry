# Assets

This directory belongs to the experiment skill `gene-imputation`.

The current refine pass upgrades the package from a static contract echo into a deterministic Tangram-shaped starter. No real biological claims should be drawn from these assets.

## What the starter now computes

- Reference cell-type centroids from raw synthetic scRNA-seq cells.
- Simplex-constrained NNLS mapping weights from shared genes in 3 spatial spots to 2 reference cell types.
- Projected expression for 6 shared genes plus 2 held-out genes, written into `imputed_expression.h5ad`.
- Machine-readable QC in `intermediate_qc.json`, including cosine similarities, mapping weights, marker recovery, and assignment entropy.

## What remains approximated

- Tangram's full PyTorch optimization is replaced with a lightweight centroid-level solver so the starter stays portable and finishes in under 5 seconds.
- Held-out validation is synthetic rather than drawn from a public paired benchmark dataset.

## What this starter does not validate

- Biological correctness on real paired scRNA-seq and spatial datasets.
- GPU-backed Tangram training, density priors, or cell-level mappings.
- Atlas selection, normalization robustness, batch correction, or probe design effects.

## Upgrade path

- Swap the synthetic inputs for a pinned paired benchmark slice while preserving the deliverable names.
- Replace the centroid solver with real `tangram-sc` preprocessing, `map_cells_to_space`, and `project_genes` when the environment can support those dependencies.
- Keep the same QC structure so held-out metrics, marker recovery, and entropy stay comparable across future versions.
