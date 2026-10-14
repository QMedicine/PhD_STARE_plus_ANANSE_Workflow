# Sample, contrast and gene-label audit

The audit covers every archived and current workflow script in this checkout, including alternative historical analyses, container/version wrappers and development tools. The inventory records each file's SHA-256, duplicate-content group and label-association review sites. This is a static/dataflow audit with targeted native tests; it does not claim that every archived alternative was rerun or that source metadata authenticates physical specimens.

## Findings and their relevance

| ID | Finding | Evidence and scope | Current workflow |
|---|---|---|---|
| L01 | The loess-offset **and TMM** ATAC preparation variants put sorted Old counts under Young headers and vice versa. | Confirmed source defect. The saved loess affinity/network comparisons establish actual historical-result impact; TMM-result impact was not separately tested. | Corrected mode uses explicit sample identities. Historical mode deliberately retains and records the mismatch. |
| L12 | The custom STARE-to-ANANSE converter averages RNA TPM across source and target states instead of selecting the state of each network. | Confirmed in `load_tpm()` and the shared `EXPR_TPM` input. This loses state-specific expression differences at edge construction; its quantitative ranking impact has not been tested. | State-specific selection is implemented and is the current default. Historical pooled-TPM outputs retain distinct provenance and can be reproduced with the explicit historical_pooled option. |
| L02 | An older STARE gene-level analysis subtracts the wrong joined columns, including nonexistent column 14. | Confirmed code defect in `STARE_work/05_STARE_ABCpp_results_analysis.sh:63–65`. | This auxiliary analysis is not called. |
| L03 | The RNA cleanup notebook can reuse gene labels from an unfiltered or different table. | Conditional identity defect: removing an interior invalid-length row shifts subsequent TPM gene labels. Running DE after the known-symbol section can instead fail because the table and label vector have different lengths. Saved execution state is unknown. | TPM aligns by gene identity; DE consumes the corresponding filtered table. These notebook sections are not executed. |
| L04 | An older peak-activity script uses version sort with lexical `join`. | Ordering-contract risk. Native GNU `join --check-order` rejects the fixture. Complete identical key lists can still work; differing key sets can lose matching rows. No saved PhD-result impact established. | Counts require identical BED4 identities and row order across samples. |
| L05 | A reverse DYNAMITE summary combines Young/Old affinity and DE with an Old/Young gABC ratio. | Confirmed contrast inconsistency in a separate archived analysis. The helper always computes Old/Young, including when called by the Young-over-Old wrapper. | DYNAMITE is not part of this STARE+ANANSE workflow. |
| L06 | The archived tfactivity auxiliary-count importer zips separate sample/file lists and assigns values by row position. | Conditional risk if auxiliary list or gene ordering differs. No effect established for the supplied raw-count input route. | Per-sample featureCounts filenames, unique gene IDs and annotation identities are checked. |
| L07 | Historical plots use JASPAR-run DE with H13CORE influence. | Previously identified mixed checkpoints; inspected directions are appropriate, rather than reversed. | Both ranking and plotting use the same directional H13CORE DE table, per the user's decision. |
| L08 | Current ATAC assembly and historical RNA recount could silently overwrite duplicate count-file identities. Other ambiguous inputs were insufficiently checked. | Latent input-validation faults; no additional swap found in the completed cohort. | New guards reject duplicate/reused primary inputs, count producers, headers and gene/peak IDs; unexpected count columns; invalid counts; reserved activity-column names; and conflicting version-stripped aliases. |
| L09 | An auxiliary plotting tool could parse composite state labels incorrectly. | `A_to_B_to_C` could mean either `A → B_to_C` or `A_to_B → C`. | Ambiguous names now require explicit source and target arguments. The production R runner already receives those separately. |
| L10 | Source metadata contains age differences across modalities and four unusual strain spellings. | Itokawa RNA/ATAC HSC rows say **20 months / 10 weeks**; in-house Hi-C descriptions say approximately 24-month-old animals. Four Hi-C strain cells say C57BL/7, /8, /9 or /10. These may be source-sheet typos; biology is not independently authenticated. | Recorded as unresolved metadata qualifications. Conditions and archived metadata were not silently rewritten. |
| L11 | Gene/TF alias conversion intentionally merges some identifiers. | Historical symbol→ID dictionaries select the last row; suffix stripping and maximum duplicate-edge probability can conflate loci or motif copies. | The historical mapping policy is preserved and explicitly distinguished from accidental sample swapping. Conflicting gene-ID→symbol aliases in the TPM port now fail. |

Exact source snippets and line references for L01–L11 are in [findings.json](../evidence/audit/label_identity_audit_v01/findings.json). L12 was added after that audit; its snippet and scope are recorded below. The evidence generator and archived receipt have not been changed. Historical files remain unchanged.

## Historical TPM state selection (L12)

This is separate from the ATAC column-assignment mistake in L01. ATAC values were assigned to the wrong sample labels; RNA values were instead **combined across both cell states**, losing state-specific expression information before scoring edges.

The archived [converter](../historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py), lines 53–60, contains:

```python
df = pd.read_csv(path, sep="\t", dtype=str)
gene_col = next((c for c in gene_id_cols if c in df.columns), df.columns[0])
df = df.set_index(gene_col, drop=False)
# take mean TPM across any sample columns
num = df.drop(columns=[gene_col], errors="ignore").apply(pd.to_numeric, errors="coerce")
tpm = num.mean(axis=1).fillna(0.0)
tpm.name = "tpm"
return tpm  # index = genes (IDs or symbols)
```

The [historical network launcher](../historical/scripts/07_network_conversion/01_Run_STARE_TF_Gene_Affinities_to_ananse_network_H13CORE_loess_offsets.sh) supplies the same four-sample `counts.tpm.tsv` as `EXPR_TPM` to both networks. There is no state or sample-selection argument. For a gene with source TPMs `[1, 1]` and target TPMs `[9, 9]`, both calls therefore receive mean TPM 5, rather than source mean 1 and target mean 9. Binding still differs between the networks, and directional DE still contributes downstream. Expression-scaling extrema can differ with the network gene/TF universe, but that does not recover the missing state selection.

The intended correction selects RNA sample columns by exact samplesheet `sample_name` and `cell_state`, averages replicates within each state, and uses that state's vector for both TF and target expression in its network. The multiplicative edge formula remains the same. See [the main-method specification](TF_GENE_AFFINITIES_BINDING_PROBABILITIES_AND_RANKING.md#41-the-custom-probability-formula). The old code and results remain historical resources. The implemented correction changes the scientific edge calculation and is distinguished from exact reproduction. Native contract tests passed all 28 cases in job 30414089; job 30414090 verified the unchanged multiplicative formula on state-specific expression vectors. The downstream-only corrected-order run is recorded separately in [the implementation and run note](STATE_SPECIFIC_TPM_AND_EXPRESSION_HEATMAP_v01.md).

## Corrected examples for archived alternatives

For L01, use the same identity list for headers and data, without separately sorting filenames:

```bash
count_files=()
for lab in "${LABELS[@]}"; do
    count_files+=("$TMP_COUNTS_DIR/${lab}.counts")
done
# Existing header is also written from LABELS.
paste <(cut -f4 "$BEDWIN") "${count_files[@]}"
```

For L02, the joined fields are:

```text
1:ID  2:Young_name  3:Young_chr  4:Young_activity  5:Young_contact
6:Young_distance  7:Young_Nenh  8:Old_name  9:Old_chr
10:Old_activity  11:Old_contact  12:Old_distance  13:Old_Nenh
```

The historical and corrected calculations are:

```awk
# Historical
dN=$7-$14; dAct=$4-$11; dCon=$5-$12; dDist=$6-$13;
# Corrected, preserving that analysis's Young-minus-Old direction
dN=$7-$13; dAct=$4-$10; dCon=$5-$11; dDist=$6-$12;
```

For L03, create labels from the actual table being analysed, after filtering, and use the same table for counts and lengths:

```r
analysis_table <- dat_agg_only_known_symbols  # choose explicitly
tpm_table <- analysis_table |> dplyr::filter(is.finite(Length), Length > 0)
counts <- as.matrix(tpm_table[, hsc_cols, drop = FALSE])
rownames(counts) <- tpm_table$Geneid
gene_labels <- tpm_table$gene_name
# DE uses analysis_table and its Geneid, not a reused TPM-label vector.
```

For L04, sort **all** join inputs by the same lexical key under `LC_ALL=C`, use `join --check-order`, enable `set -euo pipefail` and check key-set completeness. For L05, a reverse Young/Old gABC comparison requires `log2((young + pseudocount)/(old + pseudocount))`; the general implementation should accept source/target identities rather than age-specific argument names. These are reviewable correction examples, not replacements for the preserved archive or newly introduced analytical routes.

## Producer-to-consumer checks

* **RNA:** featureCounts retains its BAM identity in the header; all four actual producer BAM headings match their count-file sample identities. The archived import strips path/suffix text from each heading and joins annotation by `Geneid`; it does not assign a separate sorted sample list. Cleanup selects sample columns by name. Current recount aligns to selected sample/gene identities; TPM and DE retain those names.
* **ATAC:** all 36 archived CSV rows yield identical `Read1`, `Read2` and `Sample_ID` under the historical comma-split parser and a proper CSV parser. The four selected accessions agree with age and replicate metadata. The subsequent independent `*.counts` sort is the confirmed swap point. Current assembly rejects duplicate identities and mismatched peak rows before combining values.
* **Hi-C:** all twelve raw components' mate filenames, age-condition labels and replicate numbers agree with the archived experimental sheet. Historical High/Low merge scripts construct input and output names from the same condition/replicate pair. Current pooling is explicit in `pool_group`; standalone libraries retain their own metadata. Both rounded and fixed-grid counting use `match()` on library identities before assigning matrix columns; fixed-grid state means select columns by state metadata.
* **STARE:** historical Young columns 5/6 and Old columns 7/8 agree with the activity **headers**, which are precisely where L01 mislabels the values. Current columns derive from the validated manifest or explicit historical header order. Normalization preserves the sample column names. Contacts and activities join by `cell_state`; outputs retain `(mode, state)` tuples.
* **Networks/influence:** replicate affinity matrices align genes and TF columns by names before averaging. The archived influence launcher uses Young/Old DE for Old→Young and Old/Young DE for Young→Old. Current DE fits `target/source`, reverse output negates logFC once, and inference selects one DE file and one network per mode/state using explicit source/target keys. ANANSE subtracts source weights from target weights by edge identity. The first DE column is checked as a gene identifier before the named-index compatibility adaptation.
* **Plots/reports:** both plot styles use the exact DE file passed with the corresponding influence tuple, and derive legends from explicit target/source labels. Main and relaxed display variants share rankings and DE, with different display settings and separate output directories. The report-only route validates state membership; it necessarily trusts the external report manifest's stated table provenance.
* **Other archives/tools:** alternative STARE/DYNAMITE wrappers, Hi-C conversions/downsampling/TAD/QC, coverage-track scripts, motif mapping and archived ANANSE sources were screened for positional association sites. Numerical aggregate sorting alone is not evidence of a label swap. A sorted input list is safe when labels derive from each selected file or when paired metadata and paths travel together. Alternate methods and unresolved alias policies are documented rather than presumed interchangeable.

[study_metadata_checks.json](../evidence/audit/label_identity_audit_v01/study_metadata_checks.json) records 20 sample checks and 25 archived sample/design/DiffBind CSV checks. This establishes agreement with archived metadata, not physical specimen identity or provenance of the manually premerged FASTQs.

## Validation and boundaries

Native Slurm job **30385252**, frozen scratch snapshot `source_label_audit_v40`, completed with exit 0. All **24 tests passed without skips**, including adversarial shuffled file orders and rejected ambiguous inputs. Thirteen checkpoint files were byte-identical: both ATAC modes' count/base/assignment files, reconstructed RNA counts and six TPM/filtered tables. Both full ATAC normalized activity tables also matched byte for byte. Three actual study manifests passed existence and primary-input uniqueness checks. See [native receipts](../evidence/audit/label_identity_audit_v01/native_30385252/guarded_replay.json).

The final suite, including the auxiliary transition-parser correction, passed **all 25 tests without skips** in native Slurm job **30385299**, snapshot `source_label_audit_v41`. The existing Hi-C resolution check also passed and rejected an incorrect resolution. [Final test output](../evidence/audit/label_identity_audit_v01/native_30385299/contracts_stderr.txt). Those guard qualifications remain applicable to their original runs because they added rejection checks without altering accepted scientific calculations. The later L12 state-specific TPM correction changes the scientific edge calculation; its new networks and reports require separate qualification.

The inventory at the original audit revision includes 484 scripts (266 historical, 218 then-current), 452 unique contents and 1,486 review sites. The new state-selection and heatmap adapters are recorded and validated separately in the implementation note linked above. Inventory and findings generators make the counts refreshable. Each file was statically screened; production producer/consumer associations and the listed historical hazard sites received focused semantic review. A regex screen is not a proof of arbitrary script correctness, and this audit does not claim execution of all alternative workflows. Ongoing Hi-C jobs use their original frozen snapshots; they were not mutated or restarted for validation-only changes.
