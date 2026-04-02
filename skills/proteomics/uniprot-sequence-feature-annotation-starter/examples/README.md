# Examples

Fetch the canonical `P04637` entry and write a compact JSON summary:

```bash
python3 skills/proteomics/uniprot-sequence-feature-annotation-starter/scripts/fetch_uniprot_sequence_feature_summary.py \
  --accession P04637 \
  --out scratch/uniprot/p04637_sequence_features.json
```

The repository also keeps a checked-in canonical output snapshot at `assets/p04637_sequence_feature_summary.json`.
