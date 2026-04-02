# References

- Tangram `map_cells_to_space` API
  URL: https://tangram-sc.readthedocs.io/en/latest/classes/tangram.mapping_utils.map_cells_to_space.html
  Note: Official API surface for Tangram cell-to-space mapping, mapping modes, and training parameters.

- Tangram mapping tutorial
  URL: https://tangram-sc.readthedocs.io/en/latest/tutorial_link.html
  Note: Official notebook showing training-gene selection, mapping, and projected-gene inspection on spatial data.

- Tangram spatial mapping tutorial
  URL: https://tangram-sc.readthedocs.io/en/latest/tutorial_sq_link.html
  Note: Official walkthrough of `pp_adatas`, overlap genes, density priors, and annotation projection for mapping single-cell data onto spatial observations.

- Tangram GitHub repository
  URL: https://github.com/broadinstitute/Tangram
  Note: Reference implementation, tutorials, and issues for Tangram workflows.

- Tangram paper
  URL: https://www.nature.com/articles/s41592-021-01264-7
  Note: Primary paper for single-cell to spatial alignment, projection, and deconvolution with Tangram.

- Spatial integration benchmark note
  URL: https://www.nature.com/articles/s41592-022-01481-8
  Note: Nature Methods benchmarking note that positions Tangram strongly for transcript spatial-distribution prediction tasks.

- Scanpy `normalize_total`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html
  Note: Official preprocessing reference for library-size normalization.

- Scanpy `log1p`
  URL: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html
  Note: Official preprocessing reference for log-transformed expression values.

- Squidpy `spatial_neighbors`
  URL: https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.spatial_neighbors.html
  Note: Official spatial-neighbor graph API used here as the upgrade path beyond the lightweight 3-NN surrogate.

- CellTypist `annotate`
  URL: https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html
  Note: Official prediction entry point documenting best-match prediction and optional majority-voting review.
