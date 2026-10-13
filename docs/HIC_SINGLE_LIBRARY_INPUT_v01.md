# Regular paired-end Hi-C inputs

Regular Hi-C data require one samplesheet row per biological replicate/library, with `assay=HiC`, `paired=TRUE`, and both `fastq_1` and `fastq_2` populated. `sample_name`, `replicate`, `cell_state` and the other required metadata are explicit. BAM and cooler fields are empty for a FASTQ input. The optional `pool_group` column can be omitted entirely; blank values also mean a standalone library.

The normal route maps each mate pair independently through the pinned diffHic standard presplit compatibility mapper, fixes mates, marks duplicates within that library, and name-sorts its BAM. Counting passes single-library BAMs directly to diffHic without invoking `DIFFHIC_POOL_LIBRARIES`. Retained pairs feed both the rounded squareCounts audit and the separate fixed-grid contact route. Joint offsets are fitted across libraries, then normalized CPM is averaged within each cell state for STARE. State-level averaging is distinct from merging sample reads or BAMs.

Use `pool_group` only when explicitly combining multiple independently prepared component libraries from the same replicate and state. The High/Low MMP dataset is this special case. The validator rejects duplicate replicate rows without an explicit shared pool, conflicting pool metadata and collisions between a pool identity and a standalone library identity. `validated_inputs.json` records the resolved groups, member identities and whether a merge is required.

## Templates and practice data

- `samplesheet_HiC_single_library_template_v01.tsv`: generic unpooled Hi-C rows plus RNA/ATAC entries needed for joint workflow validation.
- `params_hic_standard_generic_v02.json`: generic mapping/counting entry settings and pinned image placeholders.
- `samplesheet_HiC_single_library_practice_v01.tsv`: six mate pairs in `/vast/projects/Immune_Aging/MMP_Hi_Low_combined_HSC_2025/raw_fastq`, plus the existing RNA/ATAC BAM inputs. It has no pooling column.
- `samplesheet_HiC_single_library_BAM_checkpoints_v01.tsv`: the same six identities pointing to the new alignment checkpoints.
- The versioned practice launchers and parameter files name separate scratch runs and qualified-result releases. Actual copies are kept in the authorized project run directory.

The six practice inputs were already merged manually by the user. The workflow treats each supplied pair as one complete library. Consequently, duplicate marking occurs across that supplied library; this differs from marking duplicates separately in High and Low fractions before pooling their BAMs. The two routes retain separate outputs and do not claim identity. Future regular libraries use this same unpooled route without any High/Low-specific labels or logic.

## Reference-header correction

The paired-prefix pilot exposed a real FASTA-header handling issue. Bowtie2 and BAM headers use the first whitespace-delimited sequence ID (for example `1`), but diffHic `cutGenome()` retained the complete Ensembl description (`1 dna:chromosome ...`). Both the index and FASTA already have the same 61 sequence IDs; this was not a UCSC/Ensembl chromosome-style mismatch.

The pipeline now uses the explicit `fragments_short_ids` operation to replace restriction-site sequence labels with the first FASTA token. Empty/duplicate resulting IDs fail. Checks require chromosome lengths and every cut-site start/end coordinate to remain unchanged. `reference_sequence_ids.json` records the full mapping. No FASTA, index or historical BAM is edited. The original `fragments` helper operation remains available for historical audits. The changed command invalidates the old fragment-task cache explicitly.

## Qualification boundaries

The bounded pilot selects the first 5,000 paired records from each real practice library, writes new derivatives to scratch, and tests mapping, mate conservation, six separate retained-pair outputs, zero merge tasks, rounded count auditing and two fixed-grid state contact inputs. Pilot rounded counts use a zero count filter to avoid applying a production depth filter to a tiny prefix. This is routing qualification, not full-library QC or a biological result.

Full-file checks independently verify gzip integrity, SHA-256 identities, mate names/counts and sequence/quality lengths. BAM checks require primary-mate conservation, valid split flags and secondary hard clipping. Only qualified BAMs proceed to full STARE/ANANSE analysis; the publication helper verifies independent-library identities and absence of merge tasks before copying final GRNs/tables. Reports are generated from qualified score tables and validated separately.
