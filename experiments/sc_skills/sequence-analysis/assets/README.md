# Assets

This directory belongs to the experiment skill `sequence-analysis`.

The current starter uses only synthetic DNA slices and local JSON tables. It is designed to exercise the portable workflow contract, not to stand in for a real biological assay.

## What the starter truly computes

- Windowed exact-or-approximate lookup over the bundled reference slices, including reverse-complement handling.
- Interval overlap between the winning reference hit and the bundled feature annotations.
- Variant-state inspection at tracked reference positions.
- A lightweight primer scan that enforces product size, GC%, Wallace Tm, GC clamp, homopolymer, and simple self-complement rules.

## What is approximated

- BLAST is represented by deterministic ungapped overlap scoring plus longest exact-seed ranking on tiny inputs.
- Primer3 is represented by local heuristic scoring rather than full thermodynamic modeling or genome-wide specificity checks.
- Feature annotation comes from bundled toy intervals instead of live Ensembl, GenBank, or UniProt fetches.

## What remains out of scope

- Real database lookups, transcript selection, and accession refresh.
- Off-target primer screening against a genome or transcriptome.
- Clinical or experimental interpretation of the tracked variant labels.

## Upgrade path

- Replace the synthetic slices with pinned public reference segments or benchmark amplicons.
- Swap the lookup scorer for BLAST+ or another documented aligner while keeping the same deliverable names.
- Replace the heuristic primer scan with Primer3 or Primer-BLAST and keep the QC tables explicit about which parameters changed.
