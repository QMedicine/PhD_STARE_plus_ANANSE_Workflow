process DIFFHIC_FRAGMENTS {
    label 'diffhic'
    cpus 2
    memory '16 GB'
    time '4h'
    input:
    path fasta
    output:
    path 'pair_param.rds', emit: fragments
    script:
    """
    diffhic_counts.R fragments "${fasta}" '${params.hic_restriction_sites[0]}' ${params.hic_overhang} pair_param.rds
    """
}
process DIFFHIC_PREPARE_PAIRS {
    tag meta.sample_name
    label 'diffhic'
    cpus 2
    memory '32 GB'
    time '24h'
    publishDir { "${params.outdir}/diffhic/pairs/${meta.sample_name}" }, mode: 'copy'
    input:
    tuple val(meta), path(bam)
    path fragments
    output:
    tuple val(meta), path("${meta.sample_name}.trimmed.h5"), emit: pairs
    path '*.qc.json', emit: qc
    path 'versions.txt', emit: versions
    script:
    """
    diffhic_counts.R prepare "${bam}" "${fragments}" '${meta.sample_name}' ${params.hic_mapq} \
      ${params.hic_chim_dist} ${params.hic_max_frag} ${params.hic_min_inward} ${params.hic_min_outward} '${meta.sample_name}'
    """
}
process DIFFHIC_SQUARE_COUNTS {
    label 'diffhic'
    cpus 4
    memory '64 GB'
    time '24h'
    publishDir "${params.outdir}/diffhic/counts", mode: 'copy'
    input:
    path pairs
    path fragments
    path samplesheet
    output:
    path 'interactions_*.rds', emit: matrices
    path '*.tsv', emit: tables
    path 'counting_qc.json', emit: qc
    path 'versions.txt', emit: versions
    script:
    def files = pairs.collect { pair -> '"'+pair+'"' }.join(' ')
    """
    diffhic_counts.R count "${fragments}" "${samplesheet}" '${params.resolutions}' '${params.diffhic_norm}' ${params.hic_count_filter} ${files}
    """
}
