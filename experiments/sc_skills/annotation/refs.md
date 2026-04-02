# References

- CELLxGENE Census quick start
  URL: https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_quick_start.html
  Note: Official quick start for `open_soma()`, `obs.read(...)`, and `get_anndata(...)`, including metadata columns such as `assay`, `cell_type`, `tissue`, `tissue_general`, `suspension_type`, and `disease`.

- CELLxGENE Census datasets table tutorial
  URL: https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_datasets.html
  Note: Official tutorial for loading `census["census_info"]["datasets"]`, inspecting `dataset_id`, `collection_name`, `dataset_h5ad_path`, and `dataset_total_cell_count`, and linking catalog rows to downloadable source files.

- `cellxgene_census.download_source_h5ad`
  URL: https://chanzuckerberg.github.io/cellxgene-census/_autosummary/cellxgene_census.download_source_h5ad.html
  Note: Official API reference showing the `dataset_id`, `to_path`, and `census_version` contract for pinning source H5AD downloads.

- Scanpy installation
  URL: https://scanpy.readthedocs.io/en/stable/installation.html
  Note: Official install page for `scanpy[leiden]` and dependency guidance used in the runtime bootstrap section.

- Scanpy `pp.normalize_total`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
  Note: Canonical library-size normalization anchor for raw UMI-like counts.

- Scanpy `pp.log1p`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html
  Note: Canonical log-transform anchor for count preprocessing.

- Scanpy `pp.pca`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.pca.html
  Note: Official PCA reference used before Harmony integration and label transfer.

- Scanpy `pp.highly_variable_genes`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
  Note: Optional HVG-selection anchor for larger references; targeted panels can skip this step.

- Scanpy `external.pp.harmony_integrate`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html
  Note: Official Harmony integration entry point and placement in the workflow after PCA and before neighbor graph construction.

- harmonypy repository
  URL: https://github.com/slowkow/harmonypy
  Note: Primary Python implementation and installation anchor for Harmony when used outside Scanpy wrappers.

- Scanpy `tl.ingest`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.ingest.html
  Note: Official reference for kNN-based label and embedding transfer from a reference AnnData object to a query.

- CellTypist `annotate`
  URL: https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html
  Note: Public design anchor for probability matrices, majority-voting refinement, and low-confidence review patterns.

- Azimuth `RunAzimuth`
  URL: https://satijalab.github.io/azimuth/reference/RunAzimuth.html
  Note: Public design anchor for hierarchical annotation levels and reference-mapping outputs.

- UTAG repository
  URL: https://github.com/ElementoLab/utag
  Note: Primary implementation anchor for UTAG, including `slide_key`, `max_dist`, `normalization_mode='l1_norm'`, and clustering-resolution guidance.

- UTAG paper
  URL: https://doi.org/10.1038/s41592-022-01657-2
  Note: Canonical paper for graph-aware tissue architecture discovery that motivates the niche-labeling step.

- Harmony paper
  URL: https://doi.org/10.1038/s41592-019-0619-0
  Note: Canonical paper for batch correction and integration across experiments.
