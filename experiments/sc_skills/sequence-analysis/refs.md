# References

- BLAST command-line applications user manual
  URL: https://www.ncbi.nlm.nih.gov/books/NBK569839/
  Note: Canonical NCBI reference for BLAST+ programs and the alignment task surface that this starter approximates locally.

- Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool.
  URL: https://pubmed.ncbi.nlm.nih.gov/2231712/
  Note: Canonical BLAST paper that motivates the lookup stage; the starter keeps only a tiny ungapped overlap surrogate.

- Primer3 manual
  URL: https://primer3.org/manual.html
  Note: Official reference for primer constraints such as product size, Tm, GC content, and primer-pair penalties.

- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3—new capabilities and interfaces.
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC3424584/
  Note: Canonical Primer3 paper; useful for understanding how the real tool exceeds the heuristic toy scan used here.

- Ensembl REST API sequence endpoint
  URL: https://rest.ensembl.org/documentation/info/sequence_region
  Note: Official sequence retrieval surface for upgrading the bundled toy references to pinned public slices.

- Ensembl REST API overlap endpoint
  URL: https://rest.ensembl.org/documentation/info/overlap_region
  Note: Official feature-overlap surface for replacing the bundled interval table with public feature fetches.

- NCBI feature table documentation
  URL: https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/
  Note: Canonical NCBI guidance on interval-style feature annotations, which informed the starter's explicit start/end feature rows.
