process RNA_FEATURECOUNTS {
    cpus 8
    memory '32 GB'
    time '4h'
    publishDir { "${params.outdir}/rna/featurecounts" }, mode: 'copy'
    tag "${meta.sample_name}"
    label 'featurecounts'
    input:
    tuple val(meta), path(bam)
    path gtf
    output:
    tuple val(meta), path("*.featureCounts.txt"), emit: counts
    path '*.summary', emit: summary
    script:
    def pairedArgs = meta.paired == 'TRUE' ? '-p --countReadPairs' : ''
    """
    featureCounts ${pairedArgs} -T ${task.cpus} -a "${gtf}" --extraAttributes gene_name,gene_biotype \
      -o "${meta.sample_name}.featureCounts.txt" "${bam}"
    """
}

process RNA_MERGE_RECOUNT {
    publishDir { "${params.outdir}/rna/recount" }, mode: 'copy'
    label 'tpm'
    input:
    path count_files, stageAs: 'counts/*'
    path selection
    output:
    path 'selected_counts.csv', emit: counts
    path 'recount_comparison.json', emit: comparison
    path 'cleanup_groups.tsv', emit: cleanup
    script:
    """
    merge_recount.py --selection "${selection}" --counts ${count_files.collect { '"'+it+'"' }.join(' ')} --output selected_counts.csv
    """
}

process RNA_TPM_FILTER {
    publishDir { "${params.outdir}/rna/tpm" }, mode: 'copy'
    label 'tpm'
    input:
    path counts
    path lengths
    path idmap
    output:
    path 'counts.tpm.tsv', emit: tpm
    path 'counts.counts_filtered.tsv', emit: counts
    path 'TFs.genes_filtered.txt', emit: tfs
    path 'versions.json', emit: versions
    script:
    """
    tpm_filter.py --counts "${counts}" --lengths "${lengths}" --idmap "${idmap}" \
      --min-count ${params.min_count} --min-count-tf ${params.min_count_tf} \
      --min-tpm ${params.min_tpm} --min-tpm-tf ${params.min_tpm_tf}
    """
}

process RNA_ANNOTATION_COUNTS {
    label 'tpm'
    publishDir { "${params.outdir}/rna/annotation_symbol_sum" }, mode: 'copy'
    input:
    path count_files, stageAs: 'counts/*'
    path samplesheet
    output:
    path 'annotation_counts.csv', emit: counts
    path 'annotation_idmap.tsv', emit: idmap
    path 'annotation_lengths.tsv', emit: lengths
    path 'annotation_groups.tsv', emit: groups
    path 'annotation_counts_provenance.json', emit: provenance
    script:
    """
    annotation_counts.py --samplesheet "${samplesheet}" --counts ${count_files.collect { count -> '"'+count+'"' }.join(' ')}
    """
}

process RNA_DESEQ2 {
    publishDir { "${params.outdir}/rna/de/${source}_to_${target}" }, mode: 'copy'
    cpus 1
    tag "${source}->${target}"
    label 'deseq2'
    input:
    tuple val(source), val(target)
    path counts
    path samplesheet
    output:
    tuple val(source), val(target), path('forward.deseq2.tsv'), emit: forward
    tuple val(target), val(source), path('reverse.deseq2.tsv'), emit: reverse
    path 'versions.txt', emit: versions
    script:
    """
    deseq_historical.R "${counts}" "${samplesheet}" '${source}' '${target}' forward.deseq2.tsv reverse.deseq2.tsv ${task.cpus} '${params.deseq_state_order.join(',')}'
    """
}
