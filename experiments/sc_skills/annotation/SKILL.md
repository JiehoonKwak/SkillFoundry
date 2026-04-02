---
name: annotation-portable-skill-ref
description: Use this experiment skill to annotate cell types for MERFISH, Xenium, CosMx, and seqFISH style single-cell-resolution spatial data with a concrete runbook for platform detection, CELLxGENE/Census or Discover reference selection, Harmony-style integration, and hierarchical label transfer; reroute spot-based assays to spatial-deconvolution and tissue-niche work to the tissue-niche-annotation skill.
---
## Purpose

Use this package when each observation is a cell and the dataset has per-cell spatial coordinates. The intended platforms are MERFISH, Xenium, CosMx, and seqFISH.

Do not use this package for Visium, Slide-seq, ST, or 10x Spatial Gene Expression. If the observations are spots, stop and route to `spatial-deconvolution`.

Do not use this package to annotate tissue niches or microenvironments. This skill ends at cell-type annotation. If you need downstream niche analysis, hand off to `tissue-niche-annotation`.

## What this package gives you

- A real-run execution playbook with concrete commands, file outputs, and decision rules.
- A deterministic offline surrogate runner in `scripts/run_annotation.py` for contract validation and fallback rehearsal.
- A runtime checker in `scripts/check_runtime.py`.
- A Census/Discover environment planner in `scripts/plan_census_env.py`.
- An offline reference ranker in `scripts/rank_reference_metadata.py`.
- A sample metadata-export schema in `examples/reference_metadata_export.tsv`.

## Real-run directory contract

Create one run directory and keep every intermediate artifact there:

```bash
RUN_DIR="$PWD/scratch/annotation_real_run"
mkdir -p "$RUN_DIR"/{qc,reference,deliverables}
```

The runbook below expects these key files:

- `00_runtime_check.json`
- `00b_census_env_plan.json`
- `01_dataset_profile.json`
- `02_census_datasets.tsv`
- `03_reference_metadata.tsv`
- `04_reference_ranked.tsv`
- `05_reference_selection.json`
- `06_query_preprocess.h5ad`
- `07_reference_preprocess.h5ad`
- `08_label_transfer_qc.tsv`
- `annotation_table.tsv`
- `marker_evidence.md`

## Step 0. Check the runtime before touching the data

Run:

```bash
python3 scripts/check_runtime.py --out "$RUN_DIR/qc/00_runtime_check.json"
```

Inspect:

- `python_version`
- `cellxgene_census_requires_python`
- `host_python_supported_for_cellxgene_census`
- `missing_modules`
- `bootstrap_commands`
- `fallback_guidance`

Root cause of the common failure on this machine:

- the local default interpreter is Python `3.13`
- current published `cellxgene-census` builds require Python `<3.13`
- therefore `pip install cellxgene-census` can fail on the host even when the package name is correct

Decision rules:

- If `host_python_supported_for_cellxgene_census` is `true`, you may use the official Census Python route in the current environment.
- If it is `false`, do not waste time retrying the same install in Python 3.13. Move immediately to Step 0b and create a supported Python 3.10-3.12 runtime or use the Discover fallback route.
- If `scanpy.external` or `harmonypy` is missing, install them before a production run. If that is impossible, use the deterministic surrogate in `scripts/run_annotation.py` only to validate workflow shape, not final biology.

## Step 0b. Choose the Census/Discover route explicitly

Run:

```bash
python3 scripts/plan_census_env.py --out "$RUN_DIR/qc/00b_census_env_plan.json"
```

Inspect:

- `recommended_route`
- `rationale`
- `routes`

Preferred route on this machine when the host stays on Python 3.13:

```bash
conda create -y -n sc-annotation-py312 python=3.12
conda activate sc-annotation-py312
pip install --upgrade pip
pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist
```

If `python3.12` exists but `conda` does not:

```bash
python3.12 -m venv .venv-census312
source .venv-census312/bin/activate
pip install --upgrade pip
pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist
```

If neither path is possible, switch to the Discover metadata-export route:

- open the CELLxGENE Discover datasets browser in a real browser
- search for `"<species> <tissue> normal"`
- filter for matching organism, tissue, disease or healthy state, and `suspension_type=cell`
- export metadata to a TSV matching `examples/reference_metadata_export.tsv`
- rank it locally with `scripts/rank_reference_metadata.py`

## Step 1. Inspect the query dataset and classify the platform

Goal: confirm that each observation is a cell, coordinates exist, and the dataset is suitable for this workflow.

Run:

```bash
python3 - <<'PY'
import json
from pathlib import Path

import anndata as ad

run_dir = Path("scratch/annotation_real_run")
query_path = Path("query.h5ad")
adata = ad.read_h5ad(query_path, backed="r")
obs_cols = list(map(str, adata.obs.columns))
var_cols = list(map(str, adata.var.columns))
coord_pairs = [
    ("x_coord", "y_coord"),
    ("x_centroid", "y_centroid"),
    ("spatial_x", "spatial_y"),
    ("x", "y"),
]
coords = next(([x, y] for x, y in coord_pairs if x in obs_cols and y in obs_cols), [])
sample_candidates = [col for col in obs_cols if any(tok in col.lower() for tok in ["sample", "batch", "patient", "donor", "fov", "roi"])]
profile = {
    "query_path": str(query_path.resolve()),
    "n_obs": int(adata.n_obs),
    "n_vars": int(adata.n_vars),
    "obs_columns": obs_cols,
    "var_columns": var_cols,
    "coordinate_keys": coords,
    "sample_or_batch_candidates": sample_candidates,
    "has_log1p_marker": "log1p" in getattr(adata, "uns", {}),
    "platform_hint": str(adata.uns.get("spatial_platform", adata.uns.get("platform", ""))),
}
run_dir.joinpath("qc/01_dataset_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
PY
```

Capture before moving on:

- `n_obs`, `n_vars`
- coordinate columns
- obvious sample, donor, patient, FOV, ROI, or batch columns
- whether counts look raw or already transformed
- platform hints in `uns`, filenames, or user metadata

Decision rules:

- Continue only if the observations are cells and you have one usable `x/y` coordinate pair.
- Stop and reroute to `spatial-deconvolution` if the platform is Visium, Slide-seq, ST, or any spot-based assay.
- Stop if there is no plausible sample/batch column and you expect multi-sample correction.

## Step 2. Search for a reference dataset

### Route A: official Census Python path in a supported Python 3.10-3.12 runtime

Run inside the supported environment:

```bash
python - <<'PY'
import cellxgene_census

with cellxgene_census.open_soma(census_version="stable") as census:
    datasets = census["census_info"]["datasets"].read().concat().to_pandas()
datasets.to_csv("scratch/annotation_real_run/reference/02_census_datasets.tsv", sep="\t", index=False)
PY
```

Then rank candidates:

```bash
python experiments/sc_skills/annotation/scripts/rank_reference_metadata.py \
  --metadata scratch/annotation_real_run/reference/02_census_datasets.tsv \
  --species human \
  --tissue lung \
  --condition normal \
  --out scratch/annotation_real_run/reference/04_reference_ranked.tsv \
  --summary-out scratch/annotation_real_run/reference/05_reference_selection.json
```

### Route B: Discover browser / export fallback

Open the Discover datasets browser and use this search pattern:

- search string: `"<species> <tissue> normal"`
- if normal is too restrictive, retry:
  - `"<species> <tissue> healthy"`
  - `"<species> <tissue> control"`
  - `"<species> <tissue> atlas"`

Review or export these columns:

- `dataset_id`
- `dataset_title`
- `collection_name`
- `species`
- `tissue`
- `disease`
- `assay`
- `suspension_type`
- `dataset_total_cell_count`
- `has_broad_labels`
- `has_subtype_labels`
- `cell_type_obs_key`
- `subtype_obs_key`
- `is_primary_data`
- `dataset_h5ad_uri`

Save the export as `03_reference_metadata.tsv`, then rank it locally:

```bash
python experiments/sc_skills/annotation/scripts/rank_reference_metadata.py \
  --metadata scratch/annotation_real_run/reference/03_reference_metadata.tsv \
  --species human \
  --tissue lung \
  --condition normal \
  --out scratch/annotation_real_run/reference/04_reference_ranked.tsv \
  --summary-out scratch/annotation_real_run/reference/05_reference_selection.json
```

Keep the highest-ranked candidate only if it has:

- matching species
- matching tissue
- healthy or normal disease state
- `suspension_type=cell`
- usable broad cell-type labels
- enough cells for a stable reference

Reject a candidate if:

- disease state is wrong for the query
- species differs
- the reference lacks cell-type labels
- the dataset is obviously too small or poorly annotated

## Step 3. Materialize and pin the selected reference

If you are in the official Census runtime:

```bash
python - <<'PY'
import json
from pathlib import Path

import cellxgene_census

selection = json.loads(Path("scratch/annotation_real_run/reference/05_reference_selection.json").read_text())
dataset_id = selection["selected_dataset_id"]
target = Path("scratch/annotation_real_run/reference/reference.h5ad")
cellxgene_census.download_source_h5ad(
    dataset_id=dataset_id,
    to_path=str(target),
    census_version="stable",
)
print(target.resolve())
PY
```

Pin in `05_reference_selection.json`:

- `dataset_id`
- `dataset_title`
- `dataset_h5ad_uri`
- chosen broad-label column
- chosen subtype column
- runtime route used: official Census or Discover fallback

If you are using the Discover fallback and cannot materialize the H5AD locally yet:

- keep the selected row and source URI
- continue only when the source file can be downloaded on another supported machine or in a supported Python 3.10-3.12 environment

## Step 4. Preprocess query and reference consistently

Representative Scanpy path:

```python
import scanpy as sc

query = sc.read_h5ad("query.h5ad")
ref = sc.read_h5ad("scratch/annotation_real_run/reference/reference.h5ad")

shared = [gene for gene in query.var_names if gene in set(ref.var_names)]
query = query[:, shared].copy()
ref = ref[:, shared].copy()

if query.X.max() > 20:
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
if ref.X.max() > 20:
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

sc.pp.highly_variable_genes(ref, n_top_genes=min(2000, ref.n_vars), subset=True)
query = query[:, ref.var_names].copy()
sc.pp.scale(ref, max_value=10)
sc.pp.scale(query, max_value=10)
sc.pp.pca(ref)
sc.pp.pca(query)
```

Write:

- `06_query_preprocess.h5ad`
- `07_reference_preprocess.h5ad`

Stop if:

- shared genes are too few
- the reference has no usable label columns after cleanup
- query and reference require incompatible identifier systems and you cannot harmonize them

## Step 5. Integrate and transfer labels

Harmony-style integration anchor:

```python
import scanpy as sc
import scanpy.external as sce

combined = ref.concatenate(query, batch_key="dataset_role", batch_categories=["reference", "query"])
sce.pp.harmony_integrate(combined, key="dataset_role", basis="X_pca", adjusted_basis="X_pca_harmony")
```

Public API name to remember: `scanpy.external.pp.harmony_integrate`.

Reference-guided label transfer anchor:

```python
sc.pp.neighbors(ref, use_rep="X_pca")
sc.tl.ingest(query, ref, obs=["broad_label", "subtype_label"], labeling_method="knn")
```

Public API name to remember: `scanpy.tl.ingest`.

Review strategy:

- broad label first
- subtype second within the selected broad class
- retain low-confidence calls explicitly instead of force-labeling them
- compare transferred labels against marker genes and tissue plausibility

Write:

- `08_label_transfer_qc.tsv`
- `annotation_table.tsv`

`annotation_table.tsv` should contain:

- `cell_id`
- `sample_id`
- `x_coord`
- `y_coord`
- `predicted_label`
- `label_source`
- `label_confidence`
- `broad_label`
- `marker_status`
- `notes`

## Step 6. Review, stop, and hand off

Write `marker_evidence.md` with:

- the selected reference and why it was chosen
- broad-to-subtype transfer logic
- low-confidence or conflicting cells
- remaining uncertainty

Stop if:

- confidence is globally poor
- the selected reference is clearly mismatched biologically
- marker evidence contradicts the transferred labels in a systematic way
- the dataset should really be treated as spot-based

Only after annotation is credible should you move to `tissue-niche-annotation`.

## Missing-tool fallback rules

- If the current Python is `3.13`, do not retry `pip install cellxgene-census` in the same interpreter. Create a Python 3.12 environment first or switch to the Discover export route.
- If `cellxgene_census` is unavailable, use the Discover/browser export path plus `scripts/rank_reference_metadata.py`.
- If `scanpy.external` or `harmonypy` is unavailable, use the deterministic surrogate in `scripts/run_annotation.py` only for workflow rehearsal.
- If you need tissue or niche labels, stop and route to `tissue-niche-annotation`. Do not bolt that logic onto this package.

## Deterministic toy validation

The package keeps an offline surrogate path so the contract can be validated reproducibly:

```bash
python3 experiments/sc_skills/annotation/scripts/run_exercise.py --outdir scratch/annotation_toy_run
python3 -m unittest experiments.sc_skills.annotation.tests.test_contract experiments.sc_skills.annotation.tests.test_annotation -v
```

The toy run proves:

- platform detection and spot-based rejection
- reference ranking and selection logic
- Harmony-style surrogate correction
- hierarchical label transfer and marker review

The toy run does not prove:

- that a real Census or Discover search was executed
- that a real source H5AD was downloaded
- that the chosen reference is biologically correct for a real study
- anything about tissue niches or microenvironments
