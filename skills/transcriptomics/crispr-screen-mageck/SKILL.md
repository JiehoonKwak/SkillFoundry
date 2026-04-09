---
name: crispr-screen-mageck
description: Use this skill to analyze pooled CRISPR screen count data using MAGeCK-style statistics — median-ratio normalization, guide-level fold changes, and gene-level scoring via robust rank aggregation (RRA). Do not use it for arrayed screens, single-cell perturbation readouts, or base-editing tiling screens.
---

## Purpose
Identify positively and negatively selected genes from pooled CRISPR knockout/activation screen count tables, using the same statistical framework as MAGeCK (Li et al., 2014).

## When to use
- You have a guide-level count table from a pooled CRISPR screen (sgRNA counts per sample).
- You need to rank genes by depletion (negative selection) or enrichment (positive selection).
- You want a self-contained analysis without installing MAGeCK or R dependencies.
- You want MAGeCK-compatible output (gene summary TSV) for downstream tools.

## When not to use
- You have single-cell readouts (Perturb-seq) — use `perturb-seq-starter` instead.
- You need the full MAGeCK-MLE model for complex experimental designs with multiple conditions.
- You need pathway-level analysis (use MAGeCK-VISA or GSEA downstream).
- You have base-editing tiling screen data requiring variant-level resolution.

## Inputs
- Tab-separated count table: columns are `sgRNA`, `gene`, then one column per sample
- Comma-separated control sample names
- Comma-separated treatment sample names

## Outputs
- JSON payload with gene-level rankings, scores, FDR, log2 fold changes, and per-guide details
- Optional MAGeCK-compatible gene summary TSV

## Requirements
- Python 3.10+
- No external dependencies (pure stdlib: csv, json, math)

## Procedure
1. Prepare a count table in MAGeCK format (TSV with sgRNA, gene, sample columns).
2. Run:
   ```bash
   python3 skills/transcriptomics/crispr-screen-mageck/scripts/run_crispr_screen.py \
     --counts path/to/counts.tsv \
     --control "ctrl_rep1,ctrl_rep2" \
     --treatment "treatment_rep1,treatment_rep2" \
     --json-out results.json \
     --tsv-out gene_summary.tsv
   ```
3. Inspect the `gene_summary` array in the JSON output, sorted by negative selection rank.
4. Check FDR < 0.05 for statistically significant hits.
5. Review `pos` scores for enriched genes (positive selection hits).

## Validation
- Script exits successfully on the toy dataset.
- Top 4 negatively selected genes are TP53, RB1, PTEN, CDKN2A (known tumor suppressors).
- Top 3 positively selected genes are MYC, KRAS, EGFR (known oncogenes).
- Neutral controls (ACTB, GAPDH, NontargetingControl) are not significant.
- All FDR values are bounded between 0 and 1.

## Failure modes and fixes
- **Ragged rows:** Ensure every row has the same number of columns as the header.
- **Sample name mismatch:** Control/treatment names must exactly match column headers.
- **All zeros for a guide:** The pseudocount (0.5) prevents division-by-zero, but guides with zero counts across all samples provide no statistical power.
- **Single replicate:** The analysis falls back to a Poisson variance assumption. Two or more replicates per condition are strongly recommended.
- **Extreme library imbalance:** If one sample has orders of magnitude more total counts, median-ratio normalization may produce extreme size factors. Check the `size_factors` in the output.

## Safety and limits
- Statistical analysis only — does not modify input files or call external services.
- The RRA implementation is a faithful approximation of MAGeCK's algorithm but may differ slightly from the official MAGeCK C implementation on edge cases.
- For production analyses with complex designs, validate results against the official MAGeCK package.

## Example
```bash
python3 skills/transcriptomics/crispr-screen-mageck/scripts/run_crispr_screen.py \
  --counts skills/transcriptomics/crispr-screen-mageck/examples/toy_counts.tsv \
  --control "ctrl_rep1,ctrl_rep2" \
  --treatment "treatment_rep1,treatment_rep2" \
  --json-out scratch/crispr/results.json \
  --tsv-out scratch/crispr/gene_summary.tsv
```

## Provenance
- MAGeCK paper: Li et al., "MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens", Genome Biology 2014. https://doi.org/10.1186/s13059-014-0554-4
- Robust Rank Aggregation: Kolde et al., "Robust rank aggregation for gene list integration and meta-analysis", Bioinformatics 2012. https://doi.org/10.1093/bioinformatics/btr709
- MAGeCK source: https://sourceforge.net/projects/mageck/
- MAGeCK documentation: https://sourceforge.net/p/mageck/wiki/Home/

## Related skills
- `perturb-seq-starter` — for single-cell CRISPR screen readouts
- `pydeseq2-differential-expression-starter` — for differential expression analysis
- `scanpy-ranked-genes-starter` — for marker gene identification
