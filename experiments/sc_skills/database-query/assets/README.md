# Assets

This experiment package keeps only tiny synthetic lookup tables and mock source payloads. They are designed to exercise the method, not to mirror full public databases.

## What the starter truly computes

- Query cleanup and lookup-key normalization.
- Entity-type scoring across gene, protein, variant, and disease patterns.
- Fan-out over bundled alias tables and mock Ensembl, MyGene.info, UniProt, ClinVar, MedGen, and OMIM style records.
- Deterministic cross-source reconciliation with machine-readable QC outputs.

## What the starter approximates

- Source ranking is a small heuristic weight table instead of the richer relevance logic used by live services.
- Mock payloads keep only a few fields needed for tiny starter-path joins and provenance.
- Disease and variant crosswalks are reduced to compact toy concept graphs.

## What remains a surrogate

- No live network access, retry handling, pagination, or rate-limit behavior.
- No full HGVS parser, transcript remapping, or clinical evidence aggregation.
- Any accession or identifier that merely preserves public syntax inside the toy payload should be treated as a mock fixture unless replaced by a pinned live snapshot.

## Upgrade path

- Replace the mock payloads with pinned JSON responses or curated local snapshots from official APIs.
- Add richer concept graphs for transcript, protein isoform, and disease-ontology reconciliation.
- Keep the existing deliverable names stable so repository-level validation continues to pass.
