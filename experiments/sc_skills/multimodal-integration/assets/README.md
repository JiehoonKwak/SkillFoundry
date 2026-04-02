# Assets

This experiment package does not ship real multimodal datasets. The runnable starter uses only the raw synthetic counts in `examples/toy_input.json`.

## What the starter computes

- RNA library-size normalization with log transform.
- Protein CLR-like normalization.
- Chromatin TF-IDF-like normalization.
- A weighted low-rank latent space chosen by `model_choice`.
- A kNN graph, batch-mixing QC, label-consistency QC, and reference-to-query label transfer.

## What is approximated

- No true totalVI or MultiVI training occurs in this starter.
- The latent space is a deterministic weighted SVD surrogate over tiny toy inputs.
- The `.h5mu` output is a MuData-inspired HDF5 container written with `h5py` so the package stays lightweight and shell-portable.

## Upgrade path

- Swap in pinned public benchmark slices or project data while keeping the same deliverable names.
- Replace the weighted SVD block with actual `scvi-tools` model training when the runtime, dependencies, and dataset size justify it.
- Promote the HDF5 surrogate to a full `muon`/`mudata` serialization only when that dependency is available and tested in the target environment.
