# phd_stare_plus_ananse_workflow

Dataset-independent reconstruction of the PhD workflow, qualified using Old and Young female C57BL/6 murine HSCs: Hi-C contacts and ATAC activities enter STARE gABC; the preserved custom converter builds expression-informed `tf_target`/`prob` networks; ANANSE contrasts these networks and ranks TF influence.

**Development draft. End-to-end exact reproduction is not yet established.** The historical scripts, current inputs and software inventories are archived. Both ATAC column-order modes are retained at the user's request; H13CORE DE tables are used for both ranking and plotting.

Read [the reconstruction plan](docs/RECONSTRUCTION_PLAN_v01.md), [the ATAC ordering evidence](docs/ATAC_ORDERING_DECISION_v01.md), and [software provenance](docs/SOFTWARE_PROVENANCE_v01.md) before running.

The [label-identity audit](docs/LABEL_IDENTITY_AUDIT_v01.md) records additional archived column/contrast hazards, study-metadata checks, current validation guards and native byte-identical checkpoint comparisons.

[TF–gene affinities, binding probabilities and ranking](docs/TF_GENE_AFFINITIES_BINDING_PROBABILITIES_AND_RANKING.md) explains the scoring equations, illustrative contribution scenarios, scaling and interpretation of network-supported versus expression-supported TFs.

[Alternative edge weights and TF influence components](docs/ALTERNATIVE_PROB_EDGE_WEIGHTS_AND_TF_INFLUENCE_COMPONENTS.md) compares native ANANSE, binding-first expression modulation, additive scoring, motif-activity regression and proposed network-focused ranking weights. These alternatives are documented for discussion and have not been implemented.

[State-specific TPM and expression heatmap](docs/STATE_SPECIFIC_TPM_AND_EXPRESSION_HEATMAP_v01.md) records the implemented TPM-selection correction, current default, explicit historical option, shared plot scales and cached-affinity rerun. Reports now include a 90-mm Arial vector heatmap and top-20 TPM/logFC table.

Current qualification and remaining limits are in [DEVELOPMENT_STATUS_v02.md](docs/DEVELOPMENT_STATUS_v02.md). Validated historical/corrected GRNs and tables are copied to the durable project directory; the updated plotting decision and separate alternate-filter report design are described in [PLOTTING_DEFAULTS_v03.md](docs/PLOTTING_DEFAULTS_v03.md). Original saved-ranking reports and both validated ATAC ordering modes have passed native qualification (12 packages, 252 PDFs); see [current status](docs/DEVELOPMENT_STATUS_v04.md).

The samplesheet is TSV with `sample_name`, `replicate`, `species`, `assay`, `paired`, `cell_type`, `condition`, `cell_state`, `fastq_1`, `fastq_2`, `bam`, and `peak`. An additional `cooler` field represents the requested cached Hi-C entry. `source2target.tsv` contains only `source` and `target`; labels must match `cell_state` exactly. Pairing is explicit, with TRUE/FALSE values.

Current staged routes:

- `preflight`: validate sample identities, pairing, assay coverage, directions and real file paths.
- `raw_prepare`: prepare RNA/ATAC FASTQ alignment/peak checkpoints and explicitly selected standard-diffHic Hi-C BAM checkpoints, using samplesheet pairing. This pilot stage does not claim a completed downstream analysis.
- `hic_diffhic_counts`: name-sorted Hi-C BAMs → restriction assignment/pruning → squareCounts and joint loess offsets, with separately labelled experimental products.
- `cooler_estimated_ice`: prepare explicitly estimated pre-blacklist corrected matrices from merged coolers and a threshold TSV, separate from historical H5 checkpoints.
- `rna_recount`: recount STAR BAMs with featureCounts 2.0.6; preserve the authenticated historical Ensembl gene selection; reproduce mean-symbol aggregation, merged-length TPM/filtering and DESeq2 1.34.0/ashr.
- `network_replay`: restart influence from authenticated state networks and directional DE tables.
- `historical_replay`: feed saved STARE affinities, TPM and directed H13CORE DE tables through the original network converters, ANANSE influence and reports.
- `atac_historical_hic_checkpoint`: rebuild ATAC/STARE in both modes with saved RNA TPM/DE and corrected H5 checkpoints.
- `bam_historical_hic_checkpoint`: rebuild ATAC activity in both ordering modes, recount RNA, use explicitly named preserved corrected Hi-C H5 checkpoints, and rerun STARE and influence. This route does not yet prove cooler-to-ICE reproduction.

With `process_fastq=true`, RNA/ATAC FASTQs enter the downstream reconstruction routes alongside BAM inputs. All four historical single-end RNA FASTQs reproduce saved sequence/count outputs and the qualified TPM/primary DE checkpoint exactly. Synthetic paired RNA/ATAC geometry tests use generic labels and qualify mate handling; they do not establish equivalence for biological paired-end data. Real ATAC FASTQ regression remains unequal. A standard diffHic Hi-C mapper candidate passed a native split-read fixture; custom _t32 equivalence remains unresolved. Separate diffHic counting passed a two-resolution synthetic test; its restriction-rounded bins are not automatically connected to STARE. Cooler-to-ICE reconstruction and alternative normalization choices remain under development. See [Hi-C mapper and Arima notes](docs/HIC_STANDARD_PRESPLIT_AND_ARIMA_v01.md). The expired scratch ICE threshold table is an identified gap. `diffhic_norm` describes the requested future diffHic route; diffHic is not used by the historical ICE-to-STARE contact route. The separate counting stage executes this setting without changing that historical route.

All computation, images, caches, logs and intermediates use `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow`. Actual launch/manifests are copied to `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs`. The v02 launcher defaults to preflight; choose a reviewed parameter JSON with `PHDSA_PARAMS_FILE`. Slurm sends END/FAIL mail to the recorded institutional address.

```bash
cd /vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
sbatch --output=logs/preflight_%j.out --error=logs/preflight_%j.err \
  /vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/launcher_Old_and_Young_HSCs_v02.sh
```

The source, work directory and Nextflow cache must be retained for `-resume`. Do not treat existing outputs or scheduler success as proof of scientific completion. Final publication will copy validated GRNs, tables, reports, plots and provenance into separate mode/direction subfolders; development results are not automatically published as final results.

Historical scripts are copied verbatim and indexed in `historical/manifest_v01.json`, with original/resolved HPC paths and SHA-256. Supporting upstream code is archived separately. `tools/verify_historical_archive.py` checks every archived copy.

Software runs inside immutable Apptainer images. The host provides Nextflow 26.04.3 and Apptainer 1.3.5. Images and scientific payloads are excluded from Git. Build recipes, observed package inventories, checksums and validation evidence are kept in the repository.

Reports embed PDF bytes inside HTML and create local browser PDF viewers and download buttons. The user-selected branching section is adapted with generic labels and explicit colour scales. Known empty filtered subgraphs produce centred TF warning panels; other plotting errors remain fatal. Main and relaxed-filter packages are separate. No rank is interpreted as experimentally verified rejuvenation or inhibition.

Workflow labels are supplied by TSV manifests; reusable code contains no HSC/Old/Young label literals. See [generic label settings](docs/GENERIC_LABELS_v01.md) and the [quick guide](docs/QUICKGUIDE_v01.md). Specific study run files and historical scripts preserve their original labels.

For new datasets, annotation_symbol_sum can derive RNA counts, symbol maps and exon lengths without a historical selected-counts file. This is a separately named method; historical_selection remains the PhD default. See [RNA_COUNT_MODES_v01.md](docs/RNA_COUNT_MODES_v01.md).

The raw Hi-C fixed-grid connector is described in [RAW_HIC_STARE_CONNECTOR_v01.md](docs/RAW_HIC_STARE_CONNECTOR_v01.md). Current execution and plot-audit boundaries are recorded in [DEVELOPMENT_STATUS_v03.md](docs/DEVELOPMENT_STATUS_v03.md). The laptop plot_audit folder contains the requested comparisons and is excluded from publication.

Regular paired-end Hi-C libraries require no merging or pooling column. See [the standard single-library route](docs/HIC_SINGLE_LIBRARY_INPUT_v01.md) for generic templates and the six-library practice dataset.

The independent-library pilot, reference-header correction and live practice job chain are recorded in [the latest status](docs/DEVELOPMENT_STATUS_v05.md).

The [14 October Slurm update](docs/SLURM_STATUS_2026-10-14.md) records completed twelve-component computation, its project-quota publication failure, and the still-running six-library mapping route.
