# Examples

Basic ontology search:

```bash
python3 skills/systems-biology/quickgo-term-search/scripts/search_quickgo_terms.py \
  --query apoptosis \
  --limit 3
```

Write the compact JSON payload to a file:

```bash
python3 skills/systems-biology/quickgo-term-search/scripts/search_quickgo_terms.py \
  --query "DNA repair" \
  --limit 2 \
  --out scratch/quickgo/dna_repair_terms.json
```
