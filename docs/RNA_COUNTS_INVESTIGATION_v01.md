# RNA counts investigation v0.1

Native Nextflow run `rna_recount_v01` (controller 30364149) completed eight tasks using real STAR BAMs, featureCounts 2.0.6, pandas 2.3.0/Python 3.13.3, and the historical DESeq2 1.34.0/R 4.1.3 image. Operational success is separate from numerical acceptance.

The selected historical matrix has 56,098 genes × 4 HSC libraries (224,392 values). The BAM recount differs from the cleaned historical CSV in 59 values across 26 genes.

Read-only three-way diagnosis establishes:

- BAM recount versus the original raw featureCounts matrix: **zero differences** across all selected gene/sample values.
- Original raw featureCounts matrix versus historical cleaned CSV: **59 differences**, identical to the mismatch set above.

Thus alignment/BAM intake and featureCounts have reproduced the durable raw counts. The discrepancy belongs to the historical count-cleanup/export stage. Examples of cleaned values are larger than their raw values, suggesting summation of duplicate gene symbols; this remains a hypothesis pending group-level verification.

Do not use the uncorrected BAM-derived RNA output as an exact-reproduction input to GRNs. Preserve the raw recount and all diagnostics. For isolated ATAC/STARE comparisons, use the authenticated historical TPM and directed H13CORE DE checkpoints explicitly.

The duplicate-group diagnostic tests exact-case and uppercase symbol grouping, compares every selected count and records explicit retained-ID/member-ID mappings only if every value matches. Do not insert unexplained gene-specific corrections or silently copy cached counts over a failed reconstruction.
