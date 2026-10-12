include { RNA_FEATURECOUNTS; RNA_MERGE_RECOUNT; RNA_TPM_FILTER; RNA_DESEQ2 } from '../modules/rna'
workflow RNA_RECONSTRUCTION {
    take:
    manifest
    samplesheet
    source2target
    main:
    rna_samples = manifest.flatMap { m -> m.samples.findAll { r -> r.assay == 'RNA' }.collect { r -> tuple(r,file(r.bam,checkIfExists:true)) } }
    RNA_FEATURECOUNTS(rna_samples,file(params.gtf,checkIfExists:true))
    RNA_MERGE_RECOUNT(RNA_FEATURECOUNTS.out.counts.map { _meta,counts -> counts }.collect(),file(params.selected_counts,checkIfExists:true))
    RNA_TPM_FILTER(RNA_MERGE_RECOUNT.out.counts,file(params.gene_lengths,checkIfExists:true),file(params.idmap,checkIfExists:true))
    contrasts = channel.fromPath(source2target).splitCsv(header:true,sep:'\t').map { r -> tuple(r.source,r.target) }
    RNA_DESEQ2(contrasts,RNA_TPM_FILTER.out.counts.first(),samplesheet)
}
