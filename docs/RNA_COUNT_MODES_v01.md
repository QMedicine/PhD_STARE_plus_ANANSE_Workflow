# RNA count modes

`rna_counts_mode=historical_selection` is the default for exact PhD reconstruction. It uses the saved selected-ID CSV, sums matching chromosome/strand alias groups according to the recovered notebook and retains the explicitly selected ID for other duplicated symbols. The original selection values are comparison evidence, not values used to fill the new counts.

`rna_counts_mode=annotation_symbol_sum` supports new datasets without a historical selected-counts file. It sums featureCounts values across all annotated IDs sharing an uppercase gene symbol, including IDs at different loci. A lexicographically chosen ID only indexes that summed symbol. This is a distinct method; it does not claim to recover the PhD representative-ID decisions.

Missing names, dot/NA values and Ensembl IDs placed in gene_name are excluded. The helper records excluded count totals, all contributing IDs, representative IDs and exact conservation of retained annotated counts. Ambiguous semicolon-separated symbol annotations fail rather than being guessed.

If `idmap` and `gene_lengths` are omitted in annotation mode, the helper derives the map from featureCounts gene_name attributes and lengths from its merged exon lengths. The existing TPM port then averages lengths per uppercase symbol, matching its documented calculation. Explicit reference maps/lengths can still be supplied. Gene-filtered counts enter DESeq2; TF-filter files are compatibility outputs and do not filter the custom STARE network converter, which uses unfiltered TPM.

Native source_v16 passed the synthetic count-conservation contract and the four-library RNA annotation run. The historical route was separately rerun and retains byte-identical selected counts, TPM, filtered counts and primary DE. Counts and DE from the two methods must remain labelled separately.
