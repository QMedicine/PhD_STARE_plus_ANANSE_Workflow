# Development checkpoint: 13 October 2026

The reusable engine derives sample, cell-type, condition and transition labels from the TSV manifests. It contains no HSC/Old/Young label literals. Historical code and study fixtures retain their real labels. Native source_v16 passed 13 contract tests, including different cell types, non-mouse contigs, reciprocal-contrast deduplication and annotation-count conservation.

Qualified results:

| Component | Evidence and boundary |
|---|---|
| RNA from all four single-end FASTQs | raw_RNA_integration_v14 (30365466), regression 30365501: saved trimmed sequence bytes and all 78,298 raw gene counts match for every sample; cleaned counts, TPM, filtered counts and primary DE are byte-identical. Resume 30365502 reused all 16 tasks. |
| Historical ATAC assignment with saved windows | atac_saved_union_v12 (30365019), regression 30365372: all four STARE matrices have identical gene IDs, 759 numeric columns and zero numeric differences; both GRNs are byte-identical to the saved PhD networks. |
| Corrected ATAC assignment | Separately labelled GRNs/influence tables are produced. Their network hashes differ from the saved results; they are not represented as historical reproduction. |
| Influence from the saved-union run | Validation 30365485: all score components and differential-edge weights match the saved results. Ranking file bytes differ; some presentation changes occur within equal-score ties. Earlier archived-network replay had raw-score differences up to 2.22e-16; that evidence remains retained. |
| Durable final GRNs/tables | Copy 30365545: 38 source products (5,618,238,964 bytes) copied and SHA-256 verified, with two generated top-20 comparison tables, under Old_and_Young_HSCs/validated_results/atac_saved_union_v12. Qualification explicitly withholds complete reports. |
| Paired RNA/ATAC geometry | Synthetic generic-state fixtures recover all 2,000 pairs and 250–450-bp fragments. ATAC uses BAMPE and produces two summits. RNA preserves both trimmed mates. These tests do not claim equivalence for biological paired-end data. |
| Annotation-derived RNA counts | rna_annotation_mode_v16 (30367493) uses annotation_symbol_sum without a selected-counts checkpoint or cached map/lengths. Native validation 30367636 confirms exact count conservation. The separately retained legacy route (30367494) remains byte-identical to the qualified PhD RNA checkpoints. |
| Estimated ICE | cooler_estimated_ice_v15 (30366634) completed six tasks using -3/+3 thresholds. Validation 30367281 confirms finite, nonnegative, symmetric matrices. These are exploratory pre-blacklist reconstructions, not authenticated historical corrected matrices. |

The complete top-20 reports remain withheld by the user's explicit strict-failure decision. The selected notebook section renders 38/40 historical panels; SATB1 (forward) and ARX (reverse) fail after its original filters. No empty panels or relaxed thresholds have been substituted. H13CORE DE tables supply ranking and plots.

Actual single-end ATAC pilot 30368308 completed; exact regression 30368309 failed on 267/206,350 windows (maximum difference six). BAM headers establish samtools 1.20/default single-thread sorting in the original alignment, while the initial pilot used 1.22.1/32 threads. An isolated 1.20 runtime is built as raw_tools_v07 (30368384); corrected pilot 30368403 and regression 30368411 are pending. See ATAC_ALIGNMENT_PROVENANCE_v01.md.

HTTPS acquisition produced the correct file size but failed ENA MD5 and gzip CRC. A separate FTP transfer timed out; its resumed complete candidate also failed MD5. Both are retained as .partial files and never enter alignment. Independent byte diagnostics found 596 differences in three regions. Job 30368276 re-downloaded those regions from the official endpoint into a separate candidate; its complete ENA MD5 matches f747dcf29f01f4c2f82847b18da2463c and gzip CRC passes. Only this authenticated complete file enters the pilots. No cause of the earlier corruption is authenticated.

Raw Hi-C still requires an accessible copy of presplit_map_t32.py and its historical runtime provenance. The referenced shared directory denies access; no replacement mapper is silently substituted. The original cooler blacklist step also remains unresolved. Estimated ICE products stay separate from the preserved corrected-H5 route.

The primary reconstructed DE table matches the historical result exactly. The requested reverse rule negates its log2 fold changes; the independently fitted historical reverse table differs slightly. Exact historical network replay uses the explicitly supplied historical directional DE files. See RNA_NUMERICAL_EQUIVALENCE_v01.md.

All compute, work, caches, images and development logs remain in the workflow scratch root. Large inputs, matrices, images and scientific products are excluded from Git; source, historical scripts, recipes, templates and compact receipts are published through the reviewed allowlist.
