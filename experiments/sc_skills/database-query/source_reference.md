# Bioinformatics Database Query Reference Decomposition

This file is decomposition guidance for the portable experiment starter. It is not an execution contract, and it intentionally removes the original Spatial Agent private data-path assumptions and hidden wrappers.

## Task families kept in scope for the starter

- Gene lookup and alias normalization
- Protein accession and name lookup
- Variant lookup with compact ClinVar- or Ensembl-style identifiers
- Disease lookup with MedGen- or OMIM-style identifiers

## Public-source shapes that motivate the starter

- Gene and cross-reference normalization: HGNC, MyGene.info, Ensembl
- Protein lookup: UniProt
- Variant lookup: ClinVar, Ensembl variation-style records
- Disease lookup: MedGen, OMIM

## Portable starter mapping

- Replace live API calls with bundled mock JSON payloads that preserve only the fields needed for tiny deterministic runs.
- Replace private parquet databases with local alias tables and concept graphs inside `examples/toy_input.json`.
- Keep the deliverable surface small:
- `query_results.tsv`
- `source_provenance.md`
- `resolution_notes.md`

## What was intentionally dropped

- Private `data_path` runtime variables
- Hidden local database wrappers
- Large external dataset dependencies
- Unrelated task families such as miRNA, virus-host interaction, and phenotype mining that would turn this into a grab-bag template

## Promotion direction

- Swap the bundled mock payloads for pinned live-source snapshots or real API adapters.
- Add stricter variant normalization and disease-ontology reconciliation once the tiny starter path is stable.
