include { DIFFHIC_FRAGMENTS; DIFFHIC_POOL_LIBRARIES; DIFFHIC_PREPARE_PAIRS; DIFFHIC_SQUARE_COUNTS; DIFFHIC_FIXED_CONTACTS } from '../modules/diffhic'
workflow DIFFHIC_COUNTS {
    take:
    manifest
    samplesheet
    raw_bams
    main:
    if (!params.fasta || !params.diffhic_container || params.hic_restriction_sites.size()!=1) { error 'diffHic counting requires --fasta, --diffhic_container and one restriction motif' }
    samples = manifest.flatMap { m -> m.samples.findAll { r -> r.assay=='HiC' && r.bam }.collect { r ->
        tuple(r,file(r.bam,checkIfExists:true))
    } }
    libraries = samples.mix(raw_bams)
    grouped = libraries.map { meta,bam -> tuple(meta.pool_group ?: meta.sample_name,meta,bam) }.groupTuple(by:0)
    prepared = grouped.map { group,metas,bams ->
        if (bams.size()>1 && !params.hic_raw_container) { error 'Hi-C pooling requires --hic_raw_container' }
        tuple(metas[0]+[sample_name:group],bams)
    }
    DIFFHIC_POOL_LIBRARIES(prepared.filter { _meta,bams -> bams.size()>1 })
    single = prepared.filter { _meta,bams -> bams.size()==1 }.map { meta,bams -> tuple(meta,bams[0]) }
    pooled = single.mix(DIFFHIC_POOL_LIBRARIES.out.bam)
    DIFFHIC_FRAGMENTS(file(params.fasta,checkIfExists:true))
    DIFFHIC_PREPARE_PAIRS(pooled,DIFFHIC_FRAGMENTS.out.fragments.first())
    DIFFHIC_SQUARE_COUNTS(DIFFHIC_PREPARE_PAIRS.out.pairs.map { _meta,pairs -> pairs }.collect(),DIFFHIC_FRAGMENTS.out.fragments.first(),samplesheet)
    contact_ch = channel.empty()
    if (params.hic_contacts_from_diffhic) {
        DIFFHIC_FIXED_CONTACTS(DIFFHIC_PREPARE_PAIRS.out.pairs.map { _meta,pairs -> pairs }.collect(),DIFFHIC_FRAGMENTS.out.fragments.first(),samplesheet,file("${projectDir}/bin/fixed_grid_helpers.R"))
        contact_ch = DIFFHIC_FIXED_CONTACTS.out.contacts.combine(manifest).flatMap { folder,m ->
            m.samples.findAll { r -> r.assay=='HiC' }.collect { r -> r.cell_state }.unique().collect { state -> tuple(state,folder.resolve(state)) }
        }
    }
    emit:
    matrices = DIFFHIC_SQUARE_COUNTS.out.matrices
    contacts = contact_ch
}
