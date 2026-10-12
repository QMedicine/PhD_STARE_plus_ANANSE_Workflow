// Raw RNA/ATAC adapters. Native qualification is required before enabling a route.
process RNA_TRIM {
    tag meta.sample_name
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    label 'raw_tools'
    cpus 8
    memory '16 GB'
    time '4h'
    input:
    tuple val(meta), path(reads)
    path adapters
    output:
    tuple val(meta), path("${meta.sample_name}.trimmed_R*.fastq.gz"), emit: reads
    path '*.cutadapt.log', emit: log
    script:
    def readFiles = reads instanceof Collection ? reads : [reads]
    def paired = meta.paired == 'TRUE'
    if (readFiles.size() != (paired ? 2 : 1)) { error 'FASTQ count and paired field differ' }
    def readArgs = readFiles.collect { read -> '"'+read+'"' }.join(' ')
    def mateArgs = paired ? "-A file:${adapters} -p '${meta.sample_name}.trimmed_R2.fastq.gz'" : ''
    """
    cutadapt -j ${task.cpus} -q ${params.rna_trim_quality} --minimum-length ${params.rna_min_read_length} \
      -a file:"${adapters}" -o '${meta.sample_name}.trimmed_R1.fastq.gz' ${mateArgs} ${readArgs} > '${meta.sample_name}.cutadapt.log' 2>&1
    """
}

process RNA_STAR {
    tag meta.sample_name
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    label 'raw_tools'
    cpus 16
    memory '30 GB'
    time '24h'
    input:
    tuple val(meta), path(reads)
    path star_index
    output:
    tuple val(meta), path("${meta.sample_name}_Aligned.sortedByCoord.out.bam"), emit: bam
    path '*Log*.out', emit: logs
    path '*SJ.out.tab', emit: junctions
    script:
    def readFiles = reads instanceof Collection ? reads : [reads]
    def readArgs = readFiles.collect { read -> '"'+read+'"' }.join(' ')
    """
    STAR --genomeDir "${star_index}" --readFilesIn ${readArgs} --readFilesCommand zcat \
      --runThreadN ${task.cpus} --outFileNamePrefix '${meta.sample_name}_' \
      --outSAMtype BAM SortedByCoordinate --twopassMode Basic --limitBAMsortRAM ${params.star_sort_ram}
    """
}

process ATAC_ALIGN_DEDUP_FILTER {
    tag meta.sample_name
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    label 'raw_tools'
    cpus 32
    memory '30 GB'
    time '24h'
    input:
    tuple val(meta), path(reads)
    path bowtie_index
    output:
    tuple val(meta), path("${meta.sample_name}.filtered.bam"), emit: bam
    path '*.dedup_metrics.txt', emit: metrics
    path '*.bowtie2.log', emit: log
    script:
    def readFiles = reads instanceof Collection ? reads : [reads]
    def paired = meta.paired == 'TRUE'
    if (readFiles.size() != (paired ? 2 : 1)) { error 'FASTQ count and paired field differ' }
    def readArgs = paired ? "-1 '${readFiles[0]}' -2 '${readFiles[1]}'" : "-U '${readFiles[0]}'"
    // Explicit excluded contig labels preserve the legacy chrM-only study setting.
    def contigs = params.atac_excluded_contigs.join(',')
    """
    mkdir -p tmp
    bowtie2 -p ${task.cpus} --very-sensitive -x "${bowtie_index}/${params.bowtie_index_name}" ${readArgs} 2> '${meta.sample_name}.bowtie2.log' \
      | samtools view -bS - | samtools sort -@ ${task.cpus} -o aligned.sorted.bam -
    java -Djava.io.tmpdir="\$PWD/tmp" -jar /opt/raw/picard/picard.jar MarkDuplicates \
      I=aligned.sorted.bam O=dedup.bam M='${meta.sample_name}.dedup_metrics.txt' REMOVE_DUPLICATES=true
    filter_atac_contigs.py --bam dedup.bam --output '${meta.sample_name}.filtered.bam' \
      --mapq ${params.atac_mapq} --exclude '${contigs}' --mode '${params.atac_contig_filter_mode}' --threads ${task.cpus}
    """
}

process ATAC_CALL_PEAKS {
    tag meta.sample_name
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    label 'raw_tools'
    cpus 2
    memory '30 GB'
    time '12h'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path(bam), path("peaks/${meta.sample_name}_summits.bed"), emit: peaks
    path 'peaks/*', emit: files
    script:
    def geometry = meta.paired == 'TRUE' ? '-f BAMPE' : '--nomodel --shift -100 --extsize 200'
    if (!params.atac_genome_size) { error '--atac_genome_size is required for peak calling' }
    """
    macs3 callpeak -t "${bam}" ${geometry} -g ${params.atac_genome_size} -q ${params.atac_peak_qvalue} \
      -n '${meta.sample_name}' --outdir peaks
    """
}

process ATAC_RAW_BLACKLIST {
    tag meta.sample_name
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    label 'raw_tools'
    cpus 2
    memory '30 GB'
    time '12h'
    input:
    tuple val(meta), path(bam)
    path blacklist
    output:
    tuple val(meta), path("${meta.sample_name}.no_blacklist.bam"), emit: bam
    script:
    """
    bedtools intersect -v -abam "${bam}" -b "${blacklist}" | samtools sort -@ ${task.cpus} -o '${meta.sample_name}.no_blacklist.bam' -
    samtools index -@ ${task.cpus} '${meta.sample_name}.no_blacklist.bam'
    """
}
