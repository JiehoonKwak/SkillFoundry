# Example

Search ChEMBL for imatinib and save a reusable JSON snapshot:

```bash
python3 skills/drug-discovery-and-cheminformatics/chembl-molecule-search/scripts/search_chembl_molecules.py \
  --query imatinib \
  --limit 3 \
  --out skills/drug-discovery-and-cheminformatics/chembl-molecule-search/assets/imatinib_search.json
```

The output includes compact ChEMBL molecule summaries with identifiers, names, types, selected properties, and representative synonyms.
