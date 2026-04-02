---
name: chembl-molecule-search
description: Use this skill to search ChEMBL molecules by keyword and return a compact compound summary. Do not use it for local chemistry calculations or unsupported bulk export.
---

## Purpose
Search ChEMBL by free-text query and return compact compound summaries suitable for quick cheminformatics triage.

## When to use
- You need a ChEMBL identifier and basic compound summary for a known drug or molecule name.
- You want a lightweight compound lookup before deeper cheminformatics or screening workflows.

## When not to use
- You need local descriptor calculation or docking.
- You need bulk export of a large ChEMBL slice.
- You need offline execution.

## Inputs
- Query string
- Optional result limit and output path

## Outputs
- JSON payload with query metadata and compact molecule summaries

## Requirements
- Python 3.10+
- Network access to `www.ebi.ac.uk`

## Procedure
1. Run `python3 skills/drug-discovery-and-cheminformatics/chembl-molecule-search/scripts/search_chembl_molecules.py --query imatinib --limit 1 --out skills/drug-discovery-and-cheminformatics/chembl-molecule-search/assets/chembl_imatinib.json`.
2. Inspect `molecules` for `chembl_id`, `preferred_name`, `molecule_type`, `max_phase`, `canonical_smiles`, and selected properties.
3. Use the returned ChEMBL IDs in downstream compound-knowledge workflows when needed.

## Validation
- Command exits successfully.
- `molecules` contains at least one ChEMBL result for a known query.
- The first result has a non-empty `chembl_id` and molecule summary fields.

## Failure modes and fixes
- No hits: refine the query or use a more canonical drug name.
- Too many hits: lower `--limit` and inspect the top results manually.
- HTTP errors: retry later and keep the query simple for smoke verification.

## Safety and limits
- Knowledge lookup only.
- This skill does not perform synthesis planning, structure optimization, or wet-lab decision making.

## Example
- `python3 skills/drug-discovery-and-cheminformatics/chembl-molecule-search/scripts/search_chembl_molecules.py --query imatinib --limit 1`

## Provenance
- ChEMBL web services docs: https://www.ebi.ac.uk/chembl/api/data/docs
- ChEMBL API base path: https://www.ebi.ac.uk/chembl/api/data/

## Related skills
- `clinicaltrials-study-search`
