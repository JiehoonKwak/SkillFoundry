# Example Data

`toy_counts.tsv` — A synthetic CRISPR screen count table with 32 guides targeting 9 genes (including non-targeting controls).

**Expected results:**

- **Depleted (negative selection):** TP53, RB1, PTEN, CDKN2A — classic tumor suppressors whose knockout promotes tumor growth (guides lost in treatment).
- **Enriched (positive selection):** MYC, KRAS, EGFR — oncogenes whose overactivation (or guide enrichment artifacts) shows increased counts.
- **Neutral:** ACTB, GAPDH, NontargetingControl — housekeeping genes and non-targeting controls should show no significant change.

This layout mirrors a typical in vivo CRISPR screen where mice are injected with a pooled sgRNA library and guides targeting essential tumor suppressors deplete over time.
