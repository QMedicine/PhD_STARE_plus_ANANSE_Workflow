# Dataset-independent labels

Workflow code reads biological labels from `sample_name`, `cell_type`, `condition` and `cell_state` in the samplesheet. Comparison identities come from exact `source` and `target` values in source2target.tsv. Output folders, plot legends and report titles use those values. Historical scripts remain unmodified evidence; study manifests retain their real labels and input paths.

`atac_modes` defaults to `["corrected_column_order"]`. Historical label swapping is opt-in and requires `historical_header_samples` and `historical_data_samples`, each containing every ATAC sample exactly once. These explicit lists reproduce the original mismatch without hard-coding any condition or cell type in the engine. The HSC run requests both modes.

`deseq_state_order` supplies factor levels explicitly; when empty, the workflow uses first appearance in the RNA samplesheet. The HSC configuration supplies its original factor order. RNA samples retain samplesheet order, matching the archived template. Reverse comparison output negates the forward log2 fold changes.

`bin/plot_branching_template.R` adapts the exact section selected by the user: unfiltered differential network for table counts, filtered network for plots, scaled component scores in the right panel. Both directions use that section and their corresponding DE table. Its archive identity and modifications are recorded in `PLOT_TEMPLATE_PROVENANCE_v02.json`. Plot filters stay unchanged; original failures are recorded in plot_status.tsv and cause the task to fail. The user explicitly chose strict failure on 12 October 2026. `make_reports=false` allows independent network/ranking qualification while a plotting-filter revision is pending; it produces no complete scientific report. Colour limits and ticks are explicit per-comparison settings because the selected notebook section reads these variables from its interactive R session.

Use `launcher_template_v03.sh` with `PHDSA_PROJECT_RUN`, `PHDSA_SOURCE_DIR`, `PHDSA_PARAMS_FILE` and an optional `PHDSA_RUN_ID`. `container_binds` belongs to the run configuration. The WEHI scheduler profile retains institute settings.

Local checks validate renamed macrophage states, both contrast directions, explicit historical column permutations, unknown-state rejection, and the selected plot section. These contract checks do not establish scientific equivalence for a new dataset.
