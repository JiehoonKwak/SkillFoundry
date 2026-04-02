# Tests

- Smoke test: `python3 skills/scientific-knowledge/europepmc-method-triage/scripts/search_europepmc.py --query "single-cell RNA-seq" --page-size 1`
- Repository-level automated coverage: `python3 -m unittest discover -s tests -p 'test_*.py'`
