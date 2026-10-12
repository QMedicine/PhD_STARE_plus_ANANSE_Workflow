process ATAC_BLACKLIST {
    tag "${meta.sample_name}"
    label 'atac_tools'
    cpus 2
    memory '30 GB'
    time '12h'
    input:
    tuple val(meta), path(bam), path(peak)
    path blacklist
    output:
    tuple val(meta), path("${meta.sample_name}.no_blacklist.bam"), path(peak), emit: bam
    script:
    """
    bedtools intersect -v -abam "${bam}" -b "${blacklist}" | samtools sort -@ ${task.cpus} -o '${meta.sample_name}.no_blacklist.bam' -
    samtools index -@ ${task.cpus} '${meta.sample_name}.no_blacklist.bam'
    """
}

process ATAC_UNION_SUMMITS {
    label 'bedtools'
    input:
    path summits, stageAs: 'summits/*'
    path chrom_sizes
    path blacklist
    output:
    path 'union_summit.named.bed', emit: windows
    script:
    """
    cat ${summits.collect { '"'+it+'"' }.join(' ')} | awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$2+1}' | sort -k1,1 -k2,2n > all_summits.1bp.bed
    bedtools slop -i all_summits.1bp.bed -g "${chrom_sizes}" -b ${params.summit_flank} | sort -k1,1 -k2,2n > windows.bed
    bedtools subtract -A -a windows.bed -b "${blacklist}" | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,"peak_"NR}' > union_summit.named.bed
    """
}

process ATAC_WINDOW_COUNTS {
    tag "${meta.sample_name}"
    label 'bedtools'
    cpus 4
    memory '24 GB'
    time '8h'
    input:
    tuple val(meta), path(bam), path(peak)
    path windows
    output:
    path "${meta.sample_name}.window_counts.tsv", emit: counts
    script:
    """
    bedtools coverage -a "${windows}" -b "${bam}" -counts > '${meta.sample_name}.window_counts.tsv'
    """
}

process ATAC_ASSEMBLE {
    tag mode
    publishDir { "${params.outdir}/${mode}/atac" }, mode: 'copy', pattern: 'column_assignment.json'
    label 'ananse'
    input:
    val mode
    path counts, stageAs: 'counts/*'
    path samplesheet
    output:
    tuple val(mode), path('union_counts.tsv'), path('union_base.tsv'), emit: counts
    path 'column_assignment.json', emit: assignment
    script:
    """
    assemble_atac.py --samplesheet "${samplesheet}" --counts ${counts.collect { '"'+it+'"' }.join(' ')} --mode '${mode}' --output union_counts.tsv --header-order '${params.historical_header_samples.join(',')}' --data-order '${params.historical_data_samples.join(',')}'
    """
}

process ATAC_NORMALIZE {
    tag mode
    label 'r_plot'
    cpus 4
    memory '24 GB'
    time '8h'
    publishDir { "${params.outdir}/${mode}/atac" }, mode: 'copy'
    input:
    tuple val(mode), path(counts), path(base)
    output:
    tuple val(mode), path('activity.tsv'), emit: activity
    path 'loessOffsets_qc.tsv', emit: qc
    path 'versions.txt', emit: versions
    script:
    """
    normalize_atac.R "${counts}" "${base}" activity.tsv
    """
}
