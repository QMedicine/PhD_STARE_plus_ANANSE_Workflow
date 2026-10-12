# Raw ATAC alignment provenance

The supplied scripts name current module versions, and the alignment commands are commented out in the single-end branch. The saved pre-blacklist BAM provides stronger evidence of the actual earlier processing.

Its stored program records identify Bowtie2 2.5.3 (--very-sensitive, 32 threads), samtools/htslib 1.20 for SAM conversion, default single-thread coordinate sorting, Picard 3.1.1 duplicate removal, and samtools 1.20 for the MAPQ≥10/unmapped/chrM substring filter. The diagnostic's final samtools view -H record is generated while reading the header; it is not an original processing step.

The initial raw-tools v06 pilot instead used samtools 1.22.1 with a 32-thread sort. It completed but exact count regression failed: 267/206,350 windows differ by at most six counts, with 33 additional pre-blacklist alignments. This discrepancy remains recorded and is not rounded away.

raw_tools_v07 adds an isolated samtools120 wrapper with the authenticated samtools/htslib 1.20 payload. Default alignment now preserves the separate BAM conversion and single-thread sort commands; the later blacklist stage retains its separately versioned samtools 1.22.1 tools. atac_samtools and atac_sort_threads make these settings explicit. Source_v17/native pilot 30368403 completed; regression 30368411 failed on 269/206,350 windows (maximum difference six), with 32 additional pre-blacklist alignments. Matching the recovered versions did not resolve the discrepancy. A separate coordinate/read-identity diagnostic compares the saved and rebuilt BAMs; no backfilled counts or successful-equivalence claim is substituted.

Evidence: evidence/audit/atac_header_diagnostic_v01/atac_header_diagnostic.json and evidence/audit/atac_fastq_regression_v01/atac_fastq_regression.json. Program records do not by themselves establish the cause of the count differences or recover unrecorded duplicate-tie traversal state.

Native read/coordinate comparison 30370418 found 130,270 differing read names. Of these, 127,818 occur in both BAMs and all have identical sequence/quality hashes; many differ in MAPQ or alternative-hit scores. This points to alignment behavior before duplicate removal rather than a count-table formatting discrepancy. It does not prove the cause by itself.

The [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) documents a per-read seed derived from read name, sequence, quality and seed option, and SAM name truncation at whitespace. A bounded hypothesis test (30370506) re-aligned 20 differing reads under six header formats. The SRA spot/length form reproduced all 20 historical alignment payloads; alternatives matched 0–7. This is subset evidence, not complete-input qualification.

The generic default atac_fastq_header_mode=preserve keeps supplied headers. An explicit sra_spot_length mode reconstructs @ACCESSION.SPOT SPOT length=READ_LENGTH, records an audit, preserves sequence/quality lines and writes a new scratch FASTQ. It rejects non-SRA-name inputs. The authenticated ENA download is never modified. Full pilot 30370523 and strict count regression 30370526 test this hypothesis across the complete library.
