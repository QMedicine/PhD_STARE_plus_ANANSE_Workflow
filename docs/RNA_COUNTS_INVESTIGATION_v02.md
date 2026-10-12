# RNA count reconstruction: qualified duplicate aggregation

Native featureCounts 2.0.6 recounts exactly matched the original raw featureCounts output for all selected genes and samples. They initially differed from the cleaned CSV in 59 of 224,392 values across 26 genes. Summing all symbol duplicates reduced this to three values but incorrectly combined two distinct loci.

The newly recovered `historical/raw_support/DESeq2_analysis_of_Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.Rmd` defines the missing rule: sum duplicates only when their chromosome and strand tokens are identical. For distinct loci, retain one representative. The historical cleaned CSV supplies the retained IDs as an explicit checkpoint; its numeric values are comparison targets, never replacement values.

Slurm job 30364334 tested this rule on all 56,098 selected genes and four libraries: zero differences across 224,392 count values, with 311 summed groups. Evidence: `evidence/audit/rna_locus_groups_v01/locus_groups_summary.json`. This establishes the cleanup rule with retained-ID checkpoint, rather than a complete de novo proof of representative selection across all original cell types.

`merge_recount.py` applies the rule to new featureCounts output, checks annotation identities, writes the group ledger and reports comparison differences. Native Nextflow integration, TPM/filtering and DESeq2 regression remain to be qualified after this change.
