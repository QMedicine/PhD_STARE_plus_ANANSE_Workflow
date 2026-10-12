include { RNA_FEATURECOUNTS; RNA_MERGE_RECOUNT; RNA_ANNOTATION_COUNTS; RNA_TPM_FILTER; RNA_DESEQ2 } from '../modules/rna'
workflow RNA_RECONSTRUCTION {
    take:
    manifest
    samplesheet
    source2target
    raw_bam
    main:
    rna_samples = manifest.flatMap { m -> m.samples.findAll { r -> r.assay == 'RNA' && r.bam }.collect { r -> tuple(r,file(r.bam,checkIfExists:true)) } }
    all_rna = rna_samples.mix(raw_bam)
    RNA_FEATURECOUNTS(all_rna,file(params.gtf,checkIfExists:true))
    if (params.rna_counts_mode == 'historical_selection') {
        RNA_MERGE_RECOUNT(RNA_FEATURECOUNTS.out.counts.map { _meta,counts -> counts }.collect(),file(params.selected_counts,checkIfExists:true))
        expression_counts = RNA_MERGE_RECOUNT.out.counts
        expression_lengths = channel.value(file(params.gene_lengths,checkIfExists:true))
        expression_map = channel.value(file(params.idmap,checkIfExists:true))
    } else {
        RNA_ANNOTATION_COUNTS(RNA_FEATURECOUNTS.out.counts.map { _meta,counts -> counts }.collect(),samplesheet)
        expression_counts = RNA_ANNOTATION_COUNTS.out.counts
        expression_lengths = params.gene_lengths ? channel.value(file(params.gene_lengths,checkIfExists:true)) : RNA_ANNOTATION_COUNTS.out.lengths.first()
        expression_map = params.idmap ? channel.value(file(params.idmap,checkIfExists:true)) : RNA_ANNOTATION_COUNTS.out.idmap.first()
    }
    RNA_TPM_FILTER(expression_counts,expression_lengths,expression_map)
    contrasts = manifest.flatMap { m -> m.primary_transitions.collect { r -> tuple(r.source,r.target) } }
    RNA_DESEQ2(contrasts,RNA_TPM_FILTER.out.counts.first(),samplesheet)
}
