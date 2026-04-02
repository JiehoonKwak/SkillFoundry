# Examples

Run the verified toy marker-ranking workflow:

```bash
slurm/envs/scanpy/bin/python skills/transcriptomics/scanpy-ranked-genes-starter/scripts/run_scanpy_ranked_genes.py \
  --input skills/transcriptomics/scanpy-ranked-genes-starter/examples/toy_counts.tsv \
  --groups skills/transcriptomics/scanpy-ranked-genes-starter/examples/toy_groups.tsv \
  --top-n 2
```
