process DIFFHIC_FRAGMENTS {
    label 'diffhic'
    publishDir "${params.outdir}/diffhic/reference", mode: 'copy', pattern: 'reference_sequence_ids.json'
    cpus 2
    memory '32 GB'
    time '4h'
    input:
    path fasta
    output:
    path 'pair_param.rds', emit: fragments
    path 'reference_sequence_ids.json', emit: reference_ids
    script:
    """
    diffhic_counts.R fragments_short_ids "${fasta}" '${params.hic_restriction_sites[0]}' ${params.hic_overhang} pair_param.rds
    """
}

process DIFFHIC_POOL_LIBRARIES {
    tag meta.sample_name
    label 'hic_raw'
    cpus 4
    memory '16 GB'
    time '8h'
    input:
    tuple val(meta), path(bams)
    output:
    tuple val(meta), path("${meta.sample_name}.pooled.name_sorted.bam"), emit: bam
    script:
    def files = bams.collect { bam -> '"'+bam+'"' }.join(' ')
    """
    samtools merge -n -@ ${task.cpus} -o '${meta.sample_name}.pooled.name_sorted.bam' ${files}
    samtools quickcheck '${meta.sample_name}.pooled.name_sorted.bam'
    """
}

process DIFFHIC_FIXED_CONTACTS {
    label 'diffhic'
    cpus 4
    memory '64 GB'
    time '24h'
    publishDir "${params.outdir}/diffhic/fixed_grid_contacts", mode: 'copy'
    input:
    path pairs
    path fragments
    path samplesheet
    path helpers
    output:
    path 'state_contacts', emit: contacts
    path '*.tsv', emit: tables
    path 'fixed_contacts_provenance.json', emit: provenance
    path 'versions.txt', emit: versions
    script:
    def files = pairs.collect { pair -> '"'+pair+'"' }.join(' ')
    """
    diffhic_fixed_contacts.R "${fragments}" "${samplesheet}" ${params.chosen_res} '${params.chromosomes.join(',')}' '${params.diffhic_norm}' "${helpers}" ${files}
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
