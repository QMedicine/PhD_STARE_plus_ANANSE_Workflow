# Quick guide

This is a development reconstruction. It has validated count cleanup, container operations and several native Nextflow stages; it is not yet an end-to-end exact-reproduction release.

1. Create a TSV samplesheet using `templates/samplesheet_template.tsv`. Include the requested twelve columns; `cooler` is an additional column for development entry. Labels must be safe identifiers and `cell_state` must consistently describe one species, cell type and condition. Use biological replicate IDs explicitly. State labels need not refer to HSCs or ageing.
2. Create source2target.tsv with only `source` and `target`, matching samplesheet cell_state labels exactly. Reverse comparisons are enabled by default.
3. Start from `templates/params_template.json`. Select a supported stage and supply its assets and container paths. The specific HSC JSON files are fixtures with real project paths; do not treat them as generic defaults.
4. Copy source to a new versioned directory under the scratch root. Preserve its checksum manifest, retained work directory and Nextflow cache. Use `templates/launcher_template_v03.sh`, setting PHDSA_PROJECT_RUN, PHDSA_SOURCE_DIR and PHDSA_PARAMS_FILE; set PHDSA_RUN_ID for a stable run name. Submit with sbatch and explicit scratch output/error paths.
5. Inspect Slurm status, task logs, trace and numerical comparison receipts. Use PHDSA_RESUME=TRUE with the same run ID for exact resume. Published files alone are insufficient for resume or completion.

Supported draft routes:

| Stage | Inputs and output boundary |
|---|---|
| preflight | Validate samples/contrasts and files; no GRN calculations |
| rna_recount | RNA BAMs → featureCounts → retained-ID cleanup → TPM/filtering → directional DE |
| network_replay | Explicit mode/state network manifest + directional DE → influence; qualified upstream provenance must be retained |
| historical_replay | Saved state affinities + TPM + directional DE → custom GRNs → influence |
| atac_historical_hic_checkpoint | Pre-blacklist ATAC BAMs + summits, authenticated corrected H5 matrices, saved RNA TPM/DE → two ATAC modes → STARE → GRNs/influence |
| bam_historical_hic_checkpoint | Recount RNA BAMs plus ATAC BAM/summit and corrected-H5 entry → downstream GRNs/influence |

`make_reports=false` explicitly stops publication at GRNs/influence for independent numerical qualification. With reports enabled, the user-selected plotting section retains its filters and fails if any top-20 plot cannot be rendered. Complete reports are currently withheld for the HSC comparisons pending filter revision.

FASTQ processing, independent cooler-to-ICE reconstruction with estimated thresholds, and alternate normalization implementations remain under development. This draft must not silently route these inputs through a different method. All scientific calculations use Apptainer; the cluster launcher loads the pinned Nextflow and Apptainer modules.
