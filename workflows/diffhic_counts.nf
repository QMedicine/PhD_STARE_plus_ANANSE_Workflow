include { DIFFHIC_FRAGMENTS; DIFFHIC_PREPARE_PAIRS; DIFFHIC_SQUARE_COUNTS } from '../modules/diffhic'
workflow DIFFHIC_COUNTS {
    take:
    manifest
    samplesheet
    main:
    if (!params.fasta || !params.diffhic_container || params.hic_restriction_sites.size()!=1) { error 'diffHic counting requires --fasta, --diffhic_container and one restriction motif' }
    samples = manifest.flatMap { m -> m.samples.findAll { r -> r.assay=='HiC' }.collect { r ->
        if (!r.bam) { error 'diffHic counting currently requires explicit name-sorted BAM checkpoints for every Hi-C sample' }
        tuple(r,file(r.bam,checkIfExists:true))
    } }
    DIFFHIC_FRAGMENTS(file(params.fasta,checkIfExists:true))
    DIFFHIC_PREPARE_PAIRS(samples,DIFFHIC_FRAGMENTS.out.fragments.first())
    DIFFHIC_SQUARE_COUNTS(DIFFHIC_PREPARE_PAIRS.out.pairs.map { _meta,pairs -> pairs }.collect(),DIFFHIC_FRAGMENTS.out.fragments.first(),samplesheet)
    emit:
    matrices = DIFFHIC_SQUARE_COUNTS.out.matrices
}
