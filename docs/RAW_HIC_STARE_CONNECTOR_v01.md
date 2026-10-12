# Raw Hi-C qualification and fixed-grid STARE connection

The user selected direct fixed-grid counting of retained 5′ read-pair positions on 13 October 2026. Restriction-site-rounded squareCounts outputs remain separate audit products; no projection of their boundaries onto a fixed grid is performed.

The standard mapper candidate remains explicitly distinct from the unrecovered custom `_t32` helper. The historical corrected-H5 route and its existing exact GRN qualification remain available independently.

## Sample identity and pooling

The study has 12 sequenced Hi-C libraries: High and Low fractions for each of six biological replicates, three per state. The optional samplesheet column `pool_group` identifies their six pools. `replicate` retains the biological replicate number; fractions are distinguished by `sample_name`. Repeated state/assay/replicate identities are accepted only for Hi-C libraries with the same explicit pool identity and consistent state/species/condition metadata. High/Low fractions are not counted as extra biological replicates.

The mapper marks duplicates separately within each sequenced library, synchronizes mate fields and name-sorts the BAM. Pooling uses name-sorted merge and retains existing duplicate flags. It does not remark duplicates across separately prepared fractions. The synthetic pooling test uses distinct QNAME prefixes and verifies that four components yield two pools with 14 retained read pairs each.

## Contact method

diffHic 1.40.0/R 4.5.1 assigns restriction fragments and excludes marked duplicates during preparePairs. Study defaults are GATC with a four-base overhang, MAPQ≥10, chim.dist=1000, max.frag=800, min.inward=1000 and min.outward=1000. The reference and alignment chromosome naming must match.

After pruning, loadChromos/loadData retrieve the retained pairs. A positive-strand 5′ coordinate is the leftmost mapped base; a negative-strand 5′ coordinate is the rightmost mapped base, using the signed reference alignment length. These coordinates become zero-based bin starts by floor(endpoint/chosen_res)*chosen_res. Counts conserve retained pairs within the configured chromosome scope. Cis contacts are emitted as a canonical upper triangle.

Joint csaw normOffsets and calculateCPM(use.offsets=TRUE, use.norm.factors=FALSE, log=FALSE) normalize the fixed-grid count matrix across biological pools. Each state contact value is the arithmetic mean of its normalized library CPM. The resulting chr*.txt.gz files enter the same STARE_GABC process used by the historical route. A fixed task-local staging name avoids collisions between the input state contact folder and STARE's state-named output folder.

This route does not execute the historical cooler blacklist/HiCExplorer ICE reconstruction. Its method prefix is `diffhic_loess_fixed_grid__`, followed by the independently recorded ATAC assignment mode. Methods text and image identities distinguish it from the preserved corrected-H5 route. The rounded count matrices, fixed counts, log offsets and contact provenance are retained separately.

Reference: [Lun and Smyth, 2015](https://doi.org/10.1186/s12859-015-0683-0), [diffHic loadData and squareCounts reference](https://bioconductor.org/packages/release/bioc/manuals/diffHic/man/diffHic.pdf).

## Entry stages and artifacts

`hic_diffhic_counts` accepts name-sorted BAMs and, with process_fastq=true, paired Hi-C FASTQs. Setting hic_contacts_from_diffhic=true also creates state contact files. `atac_diffhic_hic_contacts` combines this contact method with ATAC reconstruction and explicit RNA TPM/DE checkpoints. `bam_diffhic_hic_contacts` additionally reconstructs RNA counts/TPM/DE. The two latter stages automatically enable fixed-grid contact generation and continue through STARE, the original custom two-column network converter and ANANSE influence.

The study's raw launch manifests are templates/samplesheet_HiC_full_raw_cohort_v01.tsv and params_HiC_full_raw_cohort_v01.json. Its downstream manifests are samplesheet_HiC_raw_BAM_checkpoints_v01.tsv and params_HiC_raw_to_STARE_v01.json. Actual copies reside in the authorized project run folder. No scientific calculation or temporary file is placed in that project folder.

## Native qualification and active cohort

| Check | Evidence |
|---|---|
| Direct fixed-grid counts, signed 5′ boundary, positive finite contacts and sample identities | 30371675 / validation 30371704; six cached tasks on resume 30372096 |
| Pooling preserves biological units and counts | 30372052; two pools, 14 retained pairs per pool |
| Real STARE execution from newly generated contacts | 30372058 / validation 30372136; four affinity matrices, three expected genes and 759 finite numeric columns, with nonzero scores; six cached tasks on resume 30372277 |
| Generic identities and pooling contract tests | 30372097; 15 native tests passed |
| Full biological raw mapping | 30370773, source_v22, hic_full_raw_cohort_v22; 12 libraries, still running at this checkpoint |
| Full compressed input hashes, gzip CRC and paired QNAME/count checks | Array 30371329, running |
| Full raw BAM conservation/flags/clipping qualification | Array 30372233, gated on successful mapping and input qualification |
| Six-pool counts → STARE → both ATAC modes → directed influence | 30372275, source_v29, hic_raw_to_stare_v29; gated on all raw BAM qualifications |
| Qualified GRN/table checksum-copy to durable storage | 30372419, gated on downstream completion |

Full-cohort completion is not claimed while these jobs are active or pending. Final qualification distinguishes internal raw-data contracts from historical custom-mapper or ICE equivalence. Complete top-20 reports remain withheld under the user's strict plotting-failure rule pending review of plot_audit.
