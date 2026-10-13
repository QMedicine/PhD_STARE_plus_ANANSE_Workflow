# Development status: independent Hi-C libraries

Regular paired-end Hi-C input now has explicit generic and practice templates without `pool_group`. Each row is processed as one library; the existing singleton branch bypasses merging. Optional pools remain available only through explicit assignments. See `HIC_SINGLE_LIBRARY_INPUT_v01.md`.

## Qualified bounded practice run

The pilot used the first 5,000 mate pairs from each of the six supplied combined libraries. Driver 30381292, source_v37, completed mapping and counting with six separate retained-pair files, two state contact inputs and zero merge tasks. All 30,000 primary mate pairs were conserved. Rounded counts were retained separately. Exact resume 30381446 reused all 16 tasks. Native container contracts passed all 18 tests in 30381131; the later report-method wording guard adds one contract; all 19 tests passed natively in 30381580 (source_v39).

This qualifies routing and input conservation for paired prefixes, not full-library results. The initial pilot 30381130 failed because cutGenome retained FASTA descriptions in its sequence labels. Diagnostic jobs 30381193 and 30381232 established that BAM/index and FASTA identifiers already match. The corrected fragment operation trims descriptions and verifies unchanged chromosome lengths and cut-site coordinates. Failed logs remain under scratch. No input, index or historical BAM was changed.

## Full six-library practice execution

Frozen source_v38 and separate scratch run IDs isolate the manually combined practice dataset from the fraction-wise run. The job chain is:

| Phase | Slurm job |
|---|---|
| Full input qualification, six libraries | 30381426, array 0–5, concurrency 2 |
| Six independent FASTQ alignments | 30381427 |
| Full BAM qualification, six libraries | 30381428, array 0–5, concurrency 2 |
| Fixed-grid contacts, STARE, both ATAC modes, both ANANSE directions | 30381429 |
| GRN/table validation and checksum publication | 30381430 |
| Main and alternate-threshold report generation/publication | 30381579 |

The full input array is running; subsequent jobs are pending successful dependencies. The six original mate pairs remain read-only at the user-provided project path. All compute, caches, intermediates and logs stay in scratch. The actual samplesheets, parameters and launcher copies are in the authorized project run directory. Qualified GRNs/tables will be copied to `validated_results/hic_single_library_to_stare_v38`; reports to `validated_results/hic_single_library_reports_v38`. The report phase uses source_v39 with explicit library-layout methods wording; its earlier pending submission 30381431 was replaced before execution.

These files were merged manually before this run. Duplicate marking therefore occurs within each supplied complete library, whereas the older 12-component route marks duplicates in High and Low separately before pooling. The routes remain separately labelled; identity is not claimed.

## Existing 12-component execution

Active mapping driver 30370773 is unchanged. BAM qualification 30372233 is unchanged. The pending, not-yet-started downstream driver 30372275 was replaced by 30381432 using the verified reference-label correction from source_v38. Publisher 30372419 now depends on that replacement; report continuation 30379823 still follows the publisher. The old queued driver was cancelled only after the publisher dependency had been rebound. Its original snapshot and logs remain preserved.

No full-library completion is claimed for either live route. Drivers and standalone qualification jobs send END/FAIL email.
