# Spatial Annotation Reference Workflow

This note is the method decomposition behind the portable `annotation` experiment package.

## Scope

Use this workflow for single-cell-resolution spatial assays where each observation is one cell, such as MERFISH, Xenium, CosMx, or seqFISH.

Do not use it for spot-based assays such as Visium, Slide-seq, ST, or 10x Spatial Gene Expression. Those belong to the deconvolution workflow.

## Real execution skeleton

1. Inspect the input spatial AnnData.
   - Confirm shape, coordinate keys, value range, and candidate sample or batch columns.
   - Decide whether counts appear raw or already normalized.

2. Search for a reference atlas in CELLxGENE/Census.
   - Start with a query shaped like `{species} {tissue} normal`.
   - Read the datasets table, then rank by species match, tissue match, healthy state, label quality, and dataset size.
   - Reject references that do not have interpretable cell-type labels.
   - If the host is Python 3.13 and the Census Python package cannot be installed, switch to the Discover browser metadata-export route and rank the exported table locally.

3. Download and pin the selected reference.
   - Record `dataset_id`, census version, output path, and chosen label columns.

4. Preprocess query and reference consistently.
   - Shared-gene intersection
   - `normalize_total`, `log1p`, PCA
   - batch or sample metadata preservation

5. Integrate and transfer labels.
   - Use Harmony-style correction or a documented fallback.
   - Transfer broad labels first, then subtype labels.
   - Keep low-confidence calls explicit.

6. Stop when the data are not trustworthy.
   - no per-cell coordinates
   - spot-based platform
   - poor gene overlap
   - missing or unreliable reference labels
   - missing supported runtime for the official Census Python API without an alternate export path

7. Hand off to the tissue niche workflow only after cell-type annotation is complete.
   - This package should stop after producing confident cell-type annotations plus QC.
   - Niche or microenvironment labeling belongs in `tissue-niche-annotation`, not here.

## Public execution anchors

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="stable") as census:
    datasets = census["census_info"]["datasets"].read().concat().to_pandas()
```

```python
cellxgene_census.download_source_h5ad(
    dataset_id="YOUR_DATASET_ID",
    to_path="reference.h5ad",
    census_version="stable",
)
```

```python
import scanpy as sc
import scanpy.external as sce

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sce.pp.harmony_integrate(adata, key="batch")
```

## Portable fallback boundary

The local runner in this package keeps the execution shape but replaces live CELLxGENE search and Harmony integration with deterministic surrogates so the package remains testable offline.

## Python 3.13 compatibility note

`cellxgene-census` is currently published for Python versions `<3.13`. On hosts whose default interpreter is Python 3.13, the exact Census Python route may be blocked until you create a dedicated Python 3.10-3.12 environment. The skill should therefore document both:

1. an exact supported-environment path for the official API, and
2. a manual/export-based fallback path that still lets an agent search and rank references using CELLxGENE Discover metadata.

## Package boundary

- This package stops after cell-type annotation plus marker/confidence review.
- Tissue niche annotation is a separate downstream workflow and belongs in `experiments/sc_skills/tissue-niche-annotation`.
