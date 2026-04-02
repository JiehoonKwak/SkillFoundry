# Tests

- Smoke test: `python3 skills/scientific-knowledge/openalex-literature-search/scripts/search_openalex.py --query "single-cell RNA-seq" --per-page 1`
- Repository-level automated coverage: `python3 -m unittest discover -s tests -p 'test_*.py'`
