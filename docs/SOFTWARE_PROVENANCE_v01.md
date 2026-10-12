# Software provenance v0.1

Historical run evidence and current environment snapshots are distinct. A snapshot does not prove that every historical task used the present packages. Verbatim inventories live in `evidence/audit/`; Conda package lists and explicit locks live in `containers/`.

| Stage | Version established | Evidence |
|---|---|---|
| RNA alignment | STAR 2.7.11b | Original STAR `@PG` records in all four HSC BAM headers |
| RNA counts | featureCounts 2.0.6 | Original featureCounts command/header in durable raw counts |
| ATAC alignment | Bowtie2 2.5.3; Picard 3.1.1; OpenJDK 21.0.2 | Archived module declarations and BAM producer records; confirm historical logs where possible |
| ATAC filtering/peaks | samtools 1.22.1; bedtools 2.31.1; MACS3 3.0.1 | Archived module declarations; audit samtools 1.23.1 is the current default and is not the historical version |
| RNA aggregation/TPM/filters | pandas 2.3.0; Python 3.13.3 | Historical tfactivity software_versions.yml and native cached-container audit |
| RNA DE | R 4.1.3; DESeq2 1.34.0; ashr | Historical software_versions.yml, archived template and actual cached SIF native inventory |
| STARE | 1.0.5; modified_STARE.sh | Original run metadata; preserved SIF SHA-256 `991095738409027e345981806f3d9ffbf6e86670791c5151dc6037083b5684ee` |
| ANANSE | 0.5.1; Python 3.10.19; numpy 1.26.4; pandas 1.5.3; scikit-learn 1.7.2 | Current named environment snapshot; historical rank regression still required |
| HiCExplorer | 3.7.6; Python 3.12.11 | Current named environment snapshot; archived scripts and native SIF test |
| cooler-env | cooler 0.10.4; Python 3.10.18 | Current named environment snapshot; historical upstream provenance still under investigation |
| hicexplorer embedded cooler | cooler 0.10.4; Python 3.12.11 | Different interpreter/dependency closure from cooler-env; retain separately |
| ATAC offsets/plots | R 4.5.1; csaw 1.42.0; SummarizedExperiment 1.38.1; data.table 1.17.8; igraph 2.2.1; ggraph 2.2.2; Cairo 1.6-5 | Observed library under `/vast/projects/Immune_Aging/R_libraries/R/x86_64-pc-linux-gnu-library/4.5.1`; packaged 104-package dependency closure; numerical offsets equal the saved matrix |
| Hi-C optional diffHic analysis | diffHic 1.40.0; edgeR 4.6.3 | Current library inventory only. The supplied STARE contact route uses ICE, not diffHic loess offsets; do not treat this optional parameter as an executed step |
| Historical tfactivity driver | Nextflow 25.04.2; tfactivity a6d7fec1ebf00ddabaaa347ed9f6025626da2697 | Historical run reports and live checkout revision |
| New driver | Nextflow module 26.04.3; Apptainer 1.3.5 | Live module catalogue; native job logs. Local lint uses Nextflow 26.04.6/Java 21.0.2 and is a separate validation environment |

R packages marked Built=4.5.2 are preserved as observed, despite the library folder being named 4.5.1. Their load warnings and plotting API warnings are evidence of potential historical drift. Exact offsets were reproduced; that does not establish exact visual equivalence for every plotting dependency.

Images stay on scratch and are excluded from Git. Original containers and named environments are unchanged. New direct-SIF images record their base checksum, recipe, archived environment/package identity, test output and final SIF checksum. Container images are private development artifacts; no image redistribution is part of the repository push.
