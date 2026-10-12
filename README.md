# phd_stare_plus_ananse_workflow

Dataset-independent reconstruction of the PhD workflow, qualified using Old and Young female C57BL/6 murine HSCs: Hi-C contacts and ATAC activities enter STARE gABC; the preserved custom converter builds expression-informed `tf_target`/`prob` networks; ANANSE contrasts these networks and ranks TF influence.

**Development draft. End-to-end exact reproduction is not yet established.** The historical scripts, current inputs and software inventories are archived. Both ATAC column-order modes are retained at the user's request; H13CORE DE tables are used for both ranking and plotting.

Read [the reconstruction plan](docs/RECONSTRUCTION_PLAN_v01.md), [the ATAC ordering evidence](docs/ATAC_ORDERING_DECISION_v01.md), and [software provenance](docs/SOFTWARE_PROVENANCE_v01.md) before running.

The samplesheet is TSV with `sample_name`, `replicate`, `species`, `assay`, `paired`, `cell_type`, `condition`, `cell_state`, `fastq_1`, `fastq_2`, `bam`, and `peak`. An additional `cooler` field represents the requested cached Hi-C entry. `source2target.tsv` contains only `source` and `target`; labels must match `cell_state` exactly. Pairing is explicit, with TRUE/FALSE values.

Current staged routes:

- `preflight`: validate sample identities, pairing, assay coverage, directions and real file paths.
- `rna_recount`: recount STAR BAMs with featureCounts 2.0.6; preserve the authenticated historical Ensembl gene selection; reproduce mean-symbol aggregation, merged-length TPM/filtering and DESeq2 1.34.0/ashr.
- `network_replay`: restart influence from authenticated state networks and directional DE tables.
- `historical_replay`: feed saved STARE affinities, TPM and directed H13CORE DE tables through the original network converters, ANANSE influence and reports.
- `atac_historical_hic_checkpoint`: rebuild ATAC/STARE in both modes with saved RNA TPM/DE and corrected H5 checkpoints.
- `bam_historical_hic_checkpoint`: rebuild ATAC activity in both ordering modes, recount RNA, use explicitly named preserved corrected Hi-C H5 checkpoints, and rerun STARE and influence. This route does not yet prove cooler-to-ICE reproduction.

FASTQ entry routes, cooler-to-ICE reconstruction and alternative normalization choices remain under development. The expired scratch ICE threshold table is an identified gap. `diffhic_norm` describes the requested future diffHic route; diffHic is not used by the historical ICE-to-STARE contact route.

All computation, images, caches, logs and intermediates use `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow`. Actual launch/manifests are copied to `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs`. The v02 launcher defaults to preflight; choose a reviewed parameter JSON with `PHDSA_PARAMS_FILE`. Slurm sends END/FAIL mail to the recorded institutional address.

```bash
cd /vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
sbatch --output=logs/preflight_%j.out --error=logs/preflight_%j.err \
  /vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs/launcher_Old_and_Young_HSCs_v02.sh
```

The source, work directory and Nextflow cache must be retained for `-resume`. Do not treat existing outputs or scheduler success as proof of scientific completion. Final publication will copy validated GRNs, tables, reports, plots and provenance into separate mode/direction subfolders; development results are not automatically published as final results.

Historical scripts are copied verbatim and indexed in `historical/manifest_v01.json`, with original/resolved HPC paths and SHA-256. Supporting upstream code is archived separately. `tools/verify_historical_archive.py` checks every archived copy.

Software runs inside immutable Apptainer images. The host provides Nextflow 26.04.3 and Apptainer 1.3.5. Images and scientific payloads are excluded from Git. Build recipes, observed package inventories, checksums and validation evidence are kept in the repository.

Reports embed PDF bytes inside HTML and create local browser PDF viewers and download buttons. The user-selected branching section is adapted with generic labels and explicit colour scales. Plotting keeps strict failure for empty filtered subgraphs; complete HSC reports await filter revision. No rank is interpreted as experimentally verified rejuvenation or inhibition.

Workflow labels are supplied by TSV manifests; reusable code contains no HSC/Old/Young label literals. See [generic label settings](docs/GENERIC_LABELS_v01.md) and the [quick guide](docs/QUICKGUIDE_v01.md). Specific study run files and historical scripts preserve their original labels.
