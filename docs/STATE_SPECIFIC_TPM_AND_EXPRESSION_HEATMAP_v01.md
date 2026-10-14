# State-specific TPM correction and TF expression heatmap

## Implemented calculation

`--tpm_expression_mode state_specific` is the current default. Before the preserved STARE converter runs, [select_state_tpm.py](../bin/select_state_tpm.py) selects RNA sample columns by exact `sample_name` and `cell_state`. The converter averages only those selected columns. Source and target networks therefore have their own mean TPM vectors. No labels are inferred from column order or age-specific filenames.

The adapter rejects duplicate headers, sample names or gene IDs, absent state/sample columns, empty RNA selections, and negative or non-finite selected TPM. It writes `tpm_selection.json` beside each network, recording state, selected RNA samples and input/output checksums. This corrects input selection; replicate-affinity averaging, per-TF binding scaling, expression scaling, the multiplicative formula and native influence calculations are unchanged:

$$
P^s_{f,g}=B^s_{f,g}(0.5+0.5T^s_f)(0.5+0.5E^s_g).
$$

`--tpm_expression_mode historical_pooled` is an explicit historical-reproduction option. It is independent of the ATAC column-order mode. The archived scripts remain verbatim; previous results are not relabelled as corrected results. Existing frozen source snapshots retain their earlier behavior. When replaying an old pooled-TPM recipe through current code, specify this option explicitly.

## Permanent report addition

[prepare_tf_expression.py](../bin/prepare_tf_expression.py) produces `top20_tf_expression.tsv` in exact top-20 ranking order. It includes source/target mean TPM, their `log2(TPM+1)` values, directional logFC, adjusted p-value and influence components. Missing expression remains `NA`; it is not invented as measured zero.

[plot_tf_expression.R](../bin/plot_tf_expression.R) renders a vector PDF titled **STARE + ANANSE Results**, subtitled **TF Rank and RNA expression**. Source and target labels are supplied as arguments. It has three columns: source expression, target expression and target/source logFC. Row labels include rank and TF. Text uses Arial at 5–7 points; the page, including legends, fits within 90 mm. Cairo rounds PDF page dimensions to points; the tested width is 89.96 mm.

Set `heatmap_tpm_limits` to `[0,10]` for `log2(TPM+1)` and `heatmap_fc_limits` to `[-6,6]` for directional logFC in the JSON parameter file; these are the shared defaults. These are display limits, not filters or scoring parameters. Both contrasts use the same limits. Colours are clipped at the limits while the accompanying HTML/TSV table retains unmodified values. Negative logFC is blue, zero is white and positive logFC is red; expression is white–green. Missing values are grey.

`--heatmap_arial_font` supplies an accessible Arial regular font file when the container lacks Arial. The study run uses Arial version 5.01.2x. The file is staged as a report input and checksum-bound in the run receipt; it is embedded in the PDF and is not redistributed in Git. The reader does not need Arial installed to view the resulting PDF. No metric-compatible font is silently substituted.

The [HTML report template](../bin/build_report.py) embeds the heatmap PDF bytes alongside the existing top-10 network overview and 20 focal TF plots. Each PDF has its own zoomable browser viewing window and download button. It also renders a compact top-20 table of rank, TF, source mean TPM, target mean TPM and directional logFC. The report is self-contained; its plots do not depend on external image URLs. Main and separately labelled relaxed-display reports receive the same heatmap and table; display filters do not recalculate scores. Each report contains 22 embedded PDFs: one heatmap, one top-10 network overview and 20 focal TF plots.

## Corrected-only cached-affinity run

Run directory: `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/runs/state_specific_tpm_v42`.

Frozen execution source: `source_state_tpm_v43`. Driver job: **30414121**. Initial driver 30414102 stopped before any workflow tasks because inherited module unloading called an unavailable `conda`; resubmission from a clean module environment resolved that initialization failure. The `historical_replay` entry stage means affinity-checkpoint intake here, not execution of the historical ATAC mode. The affinity manifest contains only `corrected_column_order` and binds the four preserved matrices from `runs/atac_saved_union_v12/results/corrected_column_order/stare`. It does not use the full raw Hi-C run's products.

The run rebuilds two state-specific multiplicative networks, native ANANSE influence in both directions, and main/relaxed-display report packages. It reuses the H13CORE directional DE tables. It does not rerun ATAC counting, STARE, RNA counting, DESeq2 or Hi-C preparation, and does not run historical comparisons or alternative edge/influence formulas. The existing positive-TF-logFC eligibility gate and equal scaled-component weights remain unchanged.

Actual launcher, samplesheet, source2target, affinity manifest, DE manifest, parameters and input checksums were copied to `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/state_specific_tpm_v42_run_files`. Final publication is targeted at `validated_results/state_specific_tpm_v42` after qualification. Computation, final qualification and publication status are recorded separately; job submission alone does not establish completion.

## Validation

* Local tests verify selection under shuffled TPM columns and interleaved state rows, rejection of ambiguous/missing inputs, and forward/reverse table orientation.
* Native contract job **30414089** passed all **28 tests without skips**. The existing Hi-C-resolution guard also passed.
* Native smoke job **30414090** passed. It used the preserved converter on a small affinity fixture: the same edge scored $1/3$ in one state and $9/14$ in the other using the correct state expression, with the multiplicative formula unchanged.
* Nextflow lint passed for the main and report-only entry points.
* The smoke PDF was visually inspected; PDFKit extracted its title and all 20 ranked TF labels as selectable text. The PDF contains an embedded `ArialMT` TrueType font and a Unicode text mapping. See [compact PDF checks](../evidence/audit/state_specific_tpm_v42/heatmap_smoke_pdf_checks.json).

The first native test submissions, 30414076 and 30414077, stopped before calculation because of inherited conflicting Apptainer modules. They were diagnosed and replaced with explicit module initialization; their logs are preserved. The execution source and native test source are frozen separately from ongoing raw Hi-C jobs, which were not restarted for this correction.
