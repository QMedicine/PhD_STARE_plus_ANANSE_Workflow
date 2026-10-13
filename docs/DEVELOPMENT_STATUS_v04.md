# Development status: selected plots and reports

The selected final force-directed plotting style is now the production default. This supersedes the earlier strict-empty plotting decision; historical evidence and notebooks remain unchanged. `PLOTTING_DEFAULTS_v03.md` records the directional defaults and the separate alternate thresholds.

## Completed and qualified

- Original saved ranking replay: driver 30379708, source_v32, four HTML packages and 84 PDFs; publisher 30379822 verified embedded bytes, exact 20+1 plot counts, ranking identity and directional node values. Durable release: `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/validated_results/selected_reports_v32`.
- Both validated ATAC ordering modes from the saved-union reconstruction: driver 30380200, source_v34, eight HTML packages and 168 PDFs. The driver qualified and checksum-copied these to `validated_results/saved_union_selected_reports_v34` under the same project run directory.
- Laptop PDFKit inspection verified selectable text in all 252 PDFs and the exact warning text and cutoff in all 11 centred warning panels. Main reports contain one top-10 PDF and 20 focal PDFs, without alternate panels. Laptop navigation is `plot_audit/selected_defaults_v32/index.html` and `plot_audit/saved_union_selected_reports_v34/index.html`; scientific report packages remain excluded from Git.
- Native container contract tests: 15 passed, job 30379005. Report exact resume: job 30380201 reused all eight tasks. Nextflow lint has no errors. All 307 historical archive copies still match their archived identities.

SATB1 and ARX remain centred warning panels in both filter variants of the original saved-ranking replay. Relaxed thresholds change display only, never the ranking or component scores. Some target genes are absent from the directional filtered H13CORE DE tables and therefore retain explicit missing-value grey fill rather than an invented fold change.

## Active full-cohort Hi-C route

The 12 Hi-C library inputs have qualified. Mapping driver 30370773 remains active. The cohort comprises two age states × three biological replicate groups × two independently prepared High/Low MMP fractions. Each fraction is mapped, mate-fixed and duplicate-marked separately. High/Low components are pooled within the explicit replicate identity, retaining their duplicate flags, to form six biological pools. diffHic retains eligible 5-prime pairs; fixed-grid counts feed joint loess-offset CPM and state-level contacts for STARE. Restriction-site-rounded squareCounts remain separate audit outputs. This is the documented standard-mapper/fixed-grid candidate method; historical custom `_t32` and ICE equivalence are not claimed.

Dependencies remain: BAM validation 30372233 → downstream STARE/ANANSE 30372275 → GRN/table qualification and publication 30372419 → updated report driver 30379823. The final report driver uses frozen source_v33, creates an explicit manifest from qualified results, renders both main and alternate reports for both ATAC modes, then validates and copies them to `validated_results/raw_hic_selected_reports_v33`. Drivers and standalone qualification jobs use END/FAIL mail. Existing active mapping and inference snapshots were not changed.

Raw Hi-C mapping, full-cohort BAM qualification and downstream scientific results are not yet complete. The earlier frozen GRN/table publisher has an obsolete “pending user plot decision” report-status message; the separately versioned report qualification receipt records eventual report completion without rewriting that historical receipt.

## Retained development failures

Source_v30 rendering stopped on colon-bearing historical DE filenames. Source_v31 stages these as `directional_de.tsv`. Source_v32 additionally corrects nested plot publication paths, makes top-10 edges fully black, retains selected labels, uses white PDF backgrounds and wraps long focal titles. Failed logs and earlier plot packages remain available; neither scientific filters nor ranking scores were silently substituted.
