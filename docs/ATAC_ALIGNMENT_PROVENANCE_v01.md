# Raw ATAC alignment provenance

The supplied scripts name current module versions, and the alignment commands are commented out in the single-end branch. The saved pre-blacklist BAM provides stronger evidence of the actual earlier processing.

Its stored program records identify Bowtie2 2.5.3 (--very-sensitive, 32 threads), samtools/htslib 1.20 for SAM conversion, default single-thread coordinate sorting, Picard 3.1.1 duplicate removal, and samtools 1.20 for the MAPQ≥10/unmapped/chrM substring filter. The diagnostic's final samtools view -H record is generated while reading the header; it is not an original processing step.

The initial raw-tools v06 pilot instead used samtools 1.22.1 with a 32-thread sort. It completed but exact count regression failed: 267/206,350 windows differ by at most six counts, with 33 additional pre-blacklist alignments. This discrepancy remains recorded and is not rounded away.

raw_tools_v07 adds an isolated samtools120 wrapper with the authenticated samtools/htslib 1.20 payload. Default alignment now preserves the separate BAM conversion and single-thread sort commands; the later blacklist stage retains its separately versioned samtools 1.22.1 tools. atac_samtools and atac_sort_threads make these settings explicit. Source_v17/native pilot 30368403 and regression 30368411 test the corrected runtime. Software/version tests pass; scientific equivalence remains pending until that comparison completes.

Evidence: evidence/audit/atac_header_diagnostic_v01/atac_header_diagnostic.json and evidence/audit/atac_fastq_regression_v01/atac_fastq_regression.json. Program records do not by themselves establish the cause of the count differences or recover unrecorded duplicate-tie traversal state.
