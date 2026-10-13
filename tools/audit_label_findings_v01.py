#!/usr/bin/env python3
"""Record reviewed findings and exact archived snippets without editing the archive."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/audit/label_identity_audit_v01'
FINDINGS=[
 ('L01','CONFIRMED_SAVED_RESULT_IMPACT','ATAC condition columns are swapped by independently sorted count filenames',[
  ('historical/scripts/05_stare_inputs/03_Add_peak_width_offsets_CPM_to_union_summit_peaks_build_bed.sh',25,33),
  ('historical/scripts/05_stare_inputs/03_Add_peak_width_offsets_CPM_to_union_summit_peaks_build_bed.sh',87,95),
  ('historical/discovered/Itokawa_et_al_2022_STARE/03_Add_peak_width_TMM_CPM_to_union_summit_peaks_build_bed.sh',87,95)],
  'The same faulty assembly is present in loess-offset and TMM variants. Saved-affinity/network regression authenticates the loess historical mode, not every alternative wrapper. Current corrected mode binds file identities explicitly; historical mode deliberately records the mismatch.'),
 ('L02','CONFIRMED_ARCHIVE_CODE_DEFECT','Old STARE gene-level join uses nonexistent and shifted columns',[
  ('historical/discovered/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/05_STARE_ABCpp_results_analysis.sh',55,70)],
  'A join of two seven-column tables has thirteen columns. Old AvgActivity, AvgContact, AvgDistance and Nenh are columns 10, 11, 12 and 13, not 11, 12, 13 and 14. This auxiliary analysis is not called by the Nextflow workflow.'),
 ('L03','CONDITIONAL_ARCHIVE_IDENTITY_DEFECT','RNA notebook reuses positional gene labels after filtering or changing tables',[
  ('historical/raw_support/DESeq2_analysis_of_Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.Rmd',476,484),
  ('historical/raw_support/DESeq2_analysis_of_Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.Rmd',548,556),
  ('historical/raw_support/DESeq2_analysis_of_Itokawa_et_al_2022_RNAseq_STAR_aligned_counts_final.Rmd',601,605)],
  'Filtering an interior invalid-length gene shifts later numeric row labels when labels are taken from the unfiltered vector. After the known-symbol section overwrites gene_labels, the DE sections select dat_agg instead: different lengths cause an error; equal lengths do not establish identity. Saved notebook execution state is unknown. Current TPM and DE use keyed rows and do not run these sections.'),
 ('L04','CONDITIONAL_ARCHIVE_SORT_CONTRACT_RISK','Version-sorted peak keys supplied to lexically ordered join',[
  ('historical/discovered/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/03_Add_activity_columns_to_union_peaks.sh',26,33),
  ('historical/discovered/Itokawa_et_al_2022_nfcore_tfactivity/STARE_work/03_Add_activity_columns_to_union_peaks.sh',51,68)],
  'peak_1, peak_2, peak_10 is version order, not lexical join order. Native GNU join --check-order rejects the fixture. Identical complete key lists can still join correctly without that option; differing key sets can lose matching rows. The archive lacks pipefail and may retain partial output on errors. No effect on the saved PhD loess activity is asserted; current preparation uses checked BED4 identities.'),
 ('L05','CONFIRMED_ARCHIVE_CONTRAST_INCONSISTENCY','DYNAMITE reverse comparison retains an Old/Young gABC ratio',[
  ('historical/discovered/Itokawa_et_al_2022_STARE/Run_manual_tfactivity_DYNAMITE_ranking/10_compute_gene_gabc_ratio.py',69,78),
  ('historical/discovered/Itokawa_et_al_2022_STARE/Run_manual_tfactivity_DYNAMITE_ranking/run_manual_tfactivity_res_100000_w_5000000_m_5000000_H13_motifs_Young_over_Old_only_known_genes.sh',198,204)],
  'The affinity ratio and DE are Young/Old, but the gABC helper always calculates Old/Young under the same reverse-comparison output prefix. This separate DYNAMITE summary route is not part of STARE+ANANSE inference.'),
 ('L06','CONDITIONAL_ARCHIVE_POSITIONAL_RISK','tfactivity extra-count import relies on parallel sample/file lists and row order',[
  ('historical/supplement/.nextflow_assets/nf-core/tfactivity/modules/local/counts/combine/templates/combine.py',10,11),
  ('historical/supplement/.nextflow_assets/nf-core/tfactivity/modules/local/counts/combine/templates/combine.py',28,30)],
  'zip does not establish sample identity or equal list lengths; to_list discards gene IDs. An effect requires misordered auxiliary files/rows and is not established for the supplied raw-count route. Current reconstruction consumes per-sample featureCounts tables with keyed sample and gene identities.'),
 ('L07','KNOWN_USER_CORRECTION','Historical plotting uses JASPAR DE and H13CORE influence',[
  ('historical/scripts/09_plots/GRN_plotting_Young2Old.Rmd',32,34),
  ('historical/scripts/09_plots/GRN_plotting_Old2Young.Rmd',34,39)],
  'These are mixed analysis checkpoints, not a demonstrated sign reversal. Current ranking and both plot styles consume the same directional H13CORE DE file, as requested.'),
 ('L08','CURRENT_GUARDS_ADDED','Duplicate count files or input paths could conceal a sample association',[],
  'ATAC assembly and historical RNA recount previously allowed duplicate dictionary keys to overwrite producers. Current validation rejects duplicate count identities, duplicate headings, duplicate peak/gene IDs, extra count columns, invalid count values, reserved ATAC annotation-column identities, conflicting version-stripped gene aliases, and shared primary read/BAM/cooler paths. Shared peak files remain supported.'),
 ('L09','CURRENT_TOOL_FIXED','Plot-rendering helper ambiguously parsed composite transition names',[],
  'The auxiliary render helper previously split A_to_B_to_C at the first separator. It now requires explicit source and target for ambiguous names. The production R plot runner already receives source and target separately.'),
 ('L10','UNRESOLVED_SOURCE_METADATA','Archived ages and mouse-strain spellings need qualified reporting',[],
  'Archived Itokawa HSC metadata specifies 20-month-old and 10-week-old RNA/ATAC animals. The in-house Hi-C description specifies approximately 24-month-old animals. Four Hi-C strain cells contain C57BL/7, C57BL/8, C57BL/9 or C57BL/10. These are recorded without silently correcting source metadata or changing conditions; their biological interpretation is not authenticated.'),
 ('L11','PRESERVED_MAPPING_POLICY','Gene/TF alias conversion is many-to-one',[
  ('historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py',69,82)],
  'The converter builds symbol-to-ID dictionaries with last-row selection, strips duplicate TF suffixes and the symbol converter takes maximum duplicate-edge probability. This can conflate distinct loci or motif copies, but is a documented part of the exact historical method, not a newly established sample-label swap. No alternative alias-selection policy is introduced in this audit.'),
]
records=[]
for identity,status,title,sites,conclusion in FINDINGS:
    snippets=[]
    for path,start,end in sites:
        lines=(ROOT/path).read_text().splitlines()
        snippets.append(dict(path=path,start=start,end=end,source='\n'.join(lines[start-1:end])))
    records.append(dict(id=identity,status=status,title=title,conclusion=conclusion,snippets=snippets))
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'findings.json').write_text(json.dumps(records,indent=2)+'\n')
print('Recorded',len(records),'reviewed findings')
