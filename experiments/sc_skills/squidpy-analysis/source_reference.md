# Squidpy Spatial Analysis

Comprehensive spatial analysis tools for single-cell and spatial transcriptomics data using Squidpy.

## When to Use Squidpy

| Analysis Type | Tool |
|---------------|------|
| Build spatial graph | `sq.gr.spatial_neighbors` |
| Cell type colocalization | `sq.gr.nhood_enrichment` |
| Co-occurrence patterns | `sq.gr.co_occurrence` |
| Spatially variable genes | `sq.gr.spatial_autocorr` |
| Point pattern analysis | `sq.gr.ripley` |
| Network topology | `sq.gr.centrality_scores` |
| Cell-cell interactions | `sq.gr.ligrec` |

---

## Workflow Overview

1. **Build spatial neighbors graph** (required first step)
2. **Analyze spatial patterns** (choose based on question)
3. **Interpret results**

---

## Step 1: Build Spatial Neighbors Graph

**Tool**: `sq.gr.spatial_neighbors`

This is required before any other Squidpy analysis.

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `coord_type` | `"visium"` (hex), `"grid"` (square), or `"generic"` (any) | `"generic"` |
| `n_neighs` | Number of neighbors for generic/KNN | 6 |
| `n_rings` | Number of hex/grid rings for visium/grid | 1 |
| `delaunay` | Use Delaunay triangulation | False |
| `radius` | Radius cutoff (0 = disabled) | 0 |

### Tips

- For **Visium data**: Use `coord_type="visium"`, `n_rings=1` (6 neighbors)
- For **generic spatial data**: Use `coord_type="generic"`, `n_neighs=6-10`
- For **Delaunay triangulation**: Set `delaunay=True` for natural neighbor connections
- For **radius-based neighbors**: Set `radius` to distance threshold (e.g., 100 microns)

---

## Step 2: Choose Analysis Based on Question

### A. Are certain cell types colocalized?

**Tool**: `sq.gr.nhood_enrichment`

Tests whether cells of one type are more/less likely to neighbor cells of another type.

**Output**:
- Z-scores: positive = enriched, negative = depleted
- Heatmap showing cell type interactions

**Example questions**:
- "Are T cells enriched near tumor cells?"
- "Which cell types cluster together?"

---

### B. Do cell types co-occur spatially?

**Tool**: `sq.gr.co_occurrence`

Measures co-occurrence probability at different spatial distances.

**Parameters**:
- `interval`: Number of distance bins (default 50)
- `n_splits`: Split for computation (default 2)

**Output**:
- Co-occurrence scores per distance
- Shows spatial range of interactions

**Example questions**:
- "At what distance do T cells and tumor cells interact?"
- "How does co-occurrence change with distance?"

---

### C. Which genes are spatially variable?

**Tool**: `sq.gr.spatial_autocorr`

Identifies genes with spatial patterns using Moran's I or Geary's C statistics.

**Parameters**:
- `mode`: `"moran"` (default) or `"geary"`
- `genes`: List of genes or "highly_variable" or None (all)
- `n_perms`: Permutations for p-value (default 100)
- `n_jobs`: Parallel jobs (default 1)

**Output**:
- Moran's I: 1 = clustered, 0 = random, -1 = dispersed
- Geary's C: 0 = clustered, 1 = random
- p-values from permutation test

**Example questions**:
- "Which genes show spatial clustering?"
- "Are marker genes spatially organized?"

---

### D. What are the spatial patterns? (Point Process)

**Tool**: `sq.gr.ripley`

Characterizes point patterns using Ripley's statistics.

**Parameters**:
- `mode`: `"F"`, `"G"`, or `"L"`
  - F: Empty space function
  - G: Nearest neighbor distribution
  - L: Ripley's L (cluster detection)
- `n_simulations`: Bootstrap simulations (default 100)

**Output**:
- L > 0: Clustering at that distance
- L < 0: Dispersion/regularity
- L = 0: Complete spatial randomness

**Example questions**:
- "Are tumor cells clustered or dispersed?"
- "At what scale do clusters form?"

---

### E. What is the network structure?

**Tool**: `sq.gr.centrality_scores`

Computes graph centrality metrics per cell type.

**Parameters**:
- `mode`: `"closeness"` or `"degree"`
  - closeness: How central in the network
  - degree: Number of connections

**Output**:
- Centrality scores per cell type
- Identifies spatially central vs peripheral populations

**Example questions**:
- "Which cell types are most central?"
- "Are tumor cells at the network periphery?"

---

### F. How do cell types interact physically?

**Tool**: `sq.gr.interaction_matrix`

Computes cell-cell contact/interaction matrix.

**Parameters**:
- `normalized`: Normalize by cell counts (default True)

**Output**:
- Interaction counts/frequencies between all cell type pairs

---

### G. What ligand-receptor interactions occur?

**Tool**: `sq.gr.ligrec`

Permutation-based ligand-receptor analysis with spatial context.

**Parameters**:
- `n_perms`: Number of permutations (default 1000)
- `threshold`: Expression threshold (default 0.01)
- `corr_method`: Multiple testing correction (default "fdr_bh")
- `gene_symbols`: Column with gene symbols (optional)

**Output**:
- Significant LR pairs between cell types
- P-values from permutation test
- Mean expression values

---

## Complete Workflow Examples

### Example 1: Tumor Microenvironment Analysis

```
1. sq.gr.spatial_neighbors(coord_type="generic", n_neighs=6)
2. sq.gr.nhood_enrichment(cluster_key="cell_type")
   → Find which immune cells infiltrate tumor
3. sq.gr.co_occurrence(cluster_key="cell_type")
   → Measure distance-dependent interactions
4. sq.gr.ligrec(cluster_key="cell_type")
   → Identify signaling between tumor and immune cells
```

### Example 2: Spatial Gene Expression Analysis

```
1. sq.gr.spatial_neighbors(coord_type="visium", n_rings=1)
2. sq.gr.spatial_autocorr(mode="moran", genes="highly_variable")
   → Find spatially variable genes
3. sq.gr.ripley(cluster_key="cell_type", mode="L")
   → Characterize spatial clustering
```

### Example 3: Tissue Organization Analysis

```
1. sq.gr.spatial_neighbors(coord_type="generic", delaunay=True)
2. sq.gr.centrality_scores(cluster_key="cell_type", mode="closeness")
   → Find central cell populations
3. sq.gr.nhood_enrichment(cluster_key="cell_type")
   → Identify neighborhood preferences
```

---

## Output Files

All tools save results to the specified `save_path`:

| Tool | Outputs |
|------|---------|
| `sq.gr.spatial_neighbors` | Updated AnnData with spatial graph |
| `sq.gr.nhood_enrichment` | Heatmap PNG, Z-scores CSV |
| `sq.gr.co_occurrence` | Co-occurrence plot PNG, scores CSV |
| `sq.gr.spatial_autocorr` | Autocorrelation results CSV |
| `sq.gr.ripley` | Ripley's statistics plot PNG, CSV |
| `sq.gr.centrality_scores` | Centrality scores CSV |
| `sq.gr.interaction_matrix` | Interaction matrix PNG, CSV |
| `sq.gr.ligrec` | LR results CSV, dotplot PNG |

---

## Tips

1. **Always run `sq.gr.spatial_neighbors` first** - Other tools depend on the spatial graph.

2. **Choose coord_type carefully**:
   - Visium: Use `"visium"` for hex grid
   - Slide-seq, MERFISH: Use `"generic"`

3. **Moran's I interpretation**:
   - High positive I → Gene is spatially clustered
   - Near zero → Random distribution
   - High negative I → Gene is spatially dispersed (rare)

4. **Neighborhood enrichment interpretation**:
   - Positive z-score → Cell types colocalize more than expected
   - Negative z-score → Cell types avoid each other
   - Near zero → Random spatial distribution

5. **For large datasets**: Use `n_jobs > 1` for parallel computation in spatial_autocorr.
