include { RNA_TRIM; RNA_STAR; ATAC_ALIGN_DEDUP_FILTER; ATAC_RAW_BLACKLIST; ATAC_CALL_PEAKS } from '../modules/raw'
workflow RAW_PREPARE {
    take:
    manifest
    main:
    if (!params.raw_tools_container) { error 'FASTQ processing requires a qualified raw-tools container' }
    rna_reads = manifest.flatMap { m -> m.samples.findAll { r -> r.assay == 'RNA' && r.fastq_1 && (!params.raw_sample_names || r.sample_name in params.raw_sample_names) }.collect { r ->
        if (!params.rna_adapters || !params.star_index) { error 'RNA FASTQs require --rna_adapters and --star_index' }
        tuple(r,[r.fastq_1,r.fastq_2].findAll { read -> read }.collect { read -> file(read,checkIfExists:true) })
    } }
    atac_reads = manifest.flatMap { m -> m.samples.findAll { r -> r.assay == 'ATAC' && r.fastq_1 && (!params.raw_sample_names || r.sample_name in params.raw_sample_names) }.collect { r ->
        if (!params.bowtie_index || !params.atac_blacklist || !params.atac_genome_size) { error 'ATAC FASTQs require --bowtie_index, --atac_blacklist and --atac_genome_size' }
        tuple(r,[r.fastq_1,r.fastq_2].findAll { read -> read }.collect { read -> file(read,checkIfExists:true) })
    } }
    RNA_TRIM(rna_reads,params.rna_adapters ? channel.value(file(params.rna_adapters,checkIfExists:true)) : channel.empty())
    RNA_STAR(RNA_TRIM.out.reads,params.star_index ? channel.value(file(params.star_index,checkIfExists:true)) : channel.empty())
    ATAC_ALIGN_DEDUP_FILTER(atac_reads,params.bowtie_index ? channel.value(file(params.bowtie_index,checkIfExists:true)) : channel.empty())
    ATAC_RAW_BLACKLIST(ATAC_ALIGN_DEDUP_FILTER.out.bam,params.atac_blacklist ? channel.value(file(params.atac_blacklist,checkIfExists:true)) : channel.empty())
    ATAC_CALL_PEAKS(ATAC_RAW_BLACKLIST.out.bam)
    emit:
    rna_bam = RNA_STAR.out.bam
    atac_bam_peaks = ATAC_CALL_PEAKS.out.peaks
}
