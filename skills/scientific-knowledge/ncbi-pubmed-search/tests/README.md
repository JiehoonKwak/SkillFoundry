# Tests

- Smoke test: `python3 skills/scientific-knowledge/ncbi-pubmed-search/scripts/search_pubmed.py --term "single-cell RNA-seq" --retmax 1`
- Repository-level coverage: `python3 -m unittest discover -s tests -p 'test_*.py'`
