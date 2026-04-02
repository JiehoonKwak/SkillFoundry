# Example

Build a small cross-source brief for a common single-cell query:

```bash
python3 skills/scientific-knowledge/multi-source-literature-brief/scripts/build_literature_brief.py \
  --query "single-cell RNA-seq" \
  --limit 2 \
  --out skills/scientific-knowledge/multi-source-literature-brief/assets/single_cell_literature_brief.json
```

The output compares normalized hits from OpenAlex, Europe PMC, Crossref, and PubMed and includes simple DOI/title overlap candidates.
