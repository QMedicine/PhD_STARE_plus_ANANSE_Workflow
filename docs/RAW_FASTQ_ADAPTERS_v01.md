# Raw FASTQ adapter draft

Generic RNA trimming/STAR alignment and ATAC alignment/deduplication/filtering/MACS3 processes are drafted in modules/raw.nf. They are wired into the raw_prepare development stage; a native single-library RNA pilot is running. This stage emits alignment/peak checkpoints and does not by itself establish downstream scientific completion. All geometry comes from the paired field, with one ordered FASTQ for FALSE and two for TRUE.

Historical RNA logs establish Cutadapt 4.8 with Python 3.9.18, quality threshold 20, minimum length 20 and the six saved adapter sequences. The current shared-system interpreter has advanced to Python 3.9.21. A separately pinned Python 3.9.18 Linux image has passed native version/platform checks; raw_tools_v05.sif snapshots the observed module payloads and passes all software tests plus a real Picard duplicate-removal fixture. FASTQ scientific regressions remain required.

STAR uses the historical twopassMode Basic, coordinate-sorted BAM output and 20-GB BAM sort limit. The preserved 60-bp index corresponds to sjdbOverhang=59. A prebuilt index is an explicit reference input; the draft does not rebuild it silently.

ATAC uses Bowtie2 2.5.3 --very-sensitive, Picard 3.1.1 with OpenJDK 21.0.2, duplicate removal, MAPQ≥10 and unmapped-read exclusion. Contig exclusions are explicit parameters. The historical fixture uses the literal chrM substring filter; a distinct RNAME mode retains SAM headers and excludes exact contigs. The differences must remain labelled rather than silently modifying the original filtering method.

The post-blacklist peak caller uses MACS3 3.0.1 and q=0.05. Single-end data uses --nomodel --shift -100 --extsize 200; paired-end data uses BAMPE. Effective genome size is required explicitly for new datasets. BAM/peak checkpoint routes retain the provided peaks; they do not silently overwrite them by recalling peaks.

Raw Hi-C reproduction remains dependent on the missing custom presplit_map_t32.py mapper and its historical Python 3.8.8/runtime provenance. The user has been asked for an accessible copy. No alternate generic aligner is substituted for that method.
