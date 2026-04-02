# Example

Search for a human BRCA1 gene summary and write a JSON snapshot:

```bash
python3 skills/genomics/ncbi-gene-search/scripts/search_ncbi_gene.py \
  --symbol BRCA1 \
  --species "homo sapiens" \
  --retmax 1 \
  --out skills/genomics/ncbi-gene-search/assets/brca1_gene_summary.json
```

The output includes the query, matched Gene IDs, and a compact `genes` list derived from NCBI Gene summaries.
