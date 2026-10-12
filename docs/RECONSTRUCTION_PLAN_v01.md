# PhD workflow reconstruction v0.1

Authoritative specification: user request and corrections, 12 October 2026.
Workflow name: `phd_stare_plus_ananse_workflow`. Historical code is preserved verbatim with SHA-256 manifests; adapters are separate files.

The repository starts from an empty directory. All compute, builds, caches, development logs and intermediates belong under `/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow`. The requested launch files and validated final outputs belong under `/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs`. Historical project inputs and named environments are read-only.

## Reconstruction stages

| Stage | Historical inputs and operation | Required verification |
|---|---|---|
| RNA alignment/counts | STAR 2.7.11b; featureCounts 2.0.6; GRCm39 Ensembl 113; single-end historical reads | Recount HSC BAMs and reconcile exported/cleaned counts without inventing duplicate-selection rules |
| TPM/filtering/DE | Archived tfactivity count helpers; pandas 2.3.0/Python 3.13.3; merged GTF gene lengths; counts sum >=50, mean TPM >=0.1; DESeq2 1.34.0/R 4.1.3/ashr | Reproduce saved TPM, retained genes and both directed contrasts |
| Hi-C | Condition pooling, cooler to H5, blacklist, hicNormalize smallest, ICE with per-matrix thresholds | Recover lost scratch threshold TSV or use authenticated stored corrected H5 checkpoint explicitly; confirm cooler provenance |
| ATAC | Pre-blacklist BAM -> blacklist removal -> MACS3 single/BAMPE; unmerged 1-bp summit windows extended 250 bp each side; bedtools alignment overlaps; csaw offsets/CPM | Run `historical_column_order` and `corrected_column_order` separately; compare against saved activity and affinities |
| STARE | Preserved modified_STARE.sh; original STARE 1.0.5 SIF; 100-kb contact bins; 5-Mb gene/enhancer windows; 0.02 threshold; all_tss; -q True; -d True; H13CORE PSEM | Compare per-replicate affinities and contact/peak identities; no GRN construction through ANANSE binding/network |
| Hybrid networks | Verbatim affinity-to-network and symbol scripts | Same four-sample mean TPM for both networks; per-TF min-max; two 0.5+0.5 expression factors; max duplicate probability |
| Influence | ANANSE 0.5.1; source/target; 100,000 pre-subtraction edges; positive target-source differences; upregulated TF candidates | Compare differential edges, components, dense-rank influence and top-20 ordering; test reverse signs |
| Plots/report | Preserve R notebook graph/aesthetics; parameterize paths/TF only; H13CORE DE for both ranking and plots | Cairo vector PDFs for top20; embedded PDF bytes, offline viewer/download; inspect text selection and layout |
| Raw entry routes | Paired Hi-C, single/paired RNA and ATAC | RNA exact and paired geometry qualified; standard Hi-C mapper candidate authorized 13 October, explicit method label and split-read qualification; missing custom _t32 equivalence unresolved |

## Decisions and limitations

1. Native audit 30363816 confirms the first three cached count rows have Old values under Young headers and vice versa. This establishes the saved matrix's labels, not which matrix produced the PhD figures. The user requests both named modes, preserved snippets and comparisons to determine provenance.
2. The user selects H13CORE DE tables for both ranking and plotting, replacing JASPAR paths in the archived notebooks. Exact visual logic is preserved; original notebook data inputs remain documented.
3. Historical scripts call windows “500 bp”, but a 1-bp summit extended 250 bp per side is 501 bp. Preserve that arithmetic for replication.
4. Current R library and Conda inventory are time-stamped snapshots, not proof of historical package versions. Historical tfactivity run records establish DESeq2 1.34.0, unlike the current shared library's 1.48.2. Several current R packages report Built=4.5.2 despite the 4.5.1 library directory. Full historical equivalence needs output regression.
5. `min_tpm_tf` filters a separate TF list in tfactivity; the historical custom network converter uses the unfiltered `counts.tpm.tsv`. Do not silently apply TF filtering to the converter.
6. The default ANANSE algorithm ranks TFs with positive target/source expression changes. Inhibition recommendations are not automatically provided by its unsigned influence output.
7. The neighbouring workflow's skills supply engineering guidance; its different scientific plan, atlas policies and checkpoints are not imported into this exact-method reconstruction.

## Acceptance

Distinguish script-byte identity, numeric artifact equality, tied-rank ordering, report changes, Linux native execution, container isolation, Slurm execution and exact resume. Publish an exact-reproduction claim only after end-to-end comparisons pass. Rebuilding software or improving cleanup is not proof of unchanged numerical results.
