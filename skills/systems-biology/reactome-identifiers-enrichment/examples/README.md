# Example

Run a compact Reactome enrichment query for two DNA-damage-related genes and save the JSON summary:

```bash
python3 skills/systems-biology/reactome-identifiers-enrichment/scripts/analyze_reactome_identifiers.py \
  --identifiers BRCA1,TP53 \
  --page-size 5 \
  --out skills/systems-biology/reactome-identifiers-enrichment/assets/brca1_tp53_enrichment.json
```

The output includes pathway summaries, species coverage, Reactome resource coverage, and any identifiers Reactome did not match.
