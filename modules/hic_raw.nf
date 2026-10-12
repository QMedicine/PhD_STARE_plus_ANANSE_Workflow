process HIC_PRESPLIT_ALIGN {
    tag meta.sample_name
    label 'hic_raw'
    publishDir { "${params.outdir}/raw/${meta.assay}/${meta.sample_name}" }, mode: 'copy'
    cpus 32
    memory '64 GB'
    time '48h'
    input:
    tuple val(meta), path(reads)
    path bowtie_index
    path upstream_source
    output:
    tuple val(meta), path("${meta.sample_name}.name_sorted.bam"), emit: bam
    path '*.json', emit: provenance
    path '*.log', emit: logs
    path '*.metrics.txt', emit: metrics
    script:
    def readFiles = reads instanceof Collection ? reads : [reads]
    if (meta.paired != 'TRUE' || readFiles.size() != 2) { error 'Hi-C presplit alignment requires two mate FASTQs' }
    """
    prepare_presplit_mapper.py --source "${upstream_source}/presplit_map.py" --output task_presplit_map.py --audit mapper_compatibility.json
    python task_presplit_map.py -G "${bowtie_index}/${params.bowtie_index_name}" \
      -1 '${readFiles[0]}' -2 '${readFiles[1]}' --sig '${params.hic_ligation_signature}' \
      --cmd 'bowtie2 -p ${task.cpus}' -o presplit.bam > presplit.log 2>&1
    mkdir -p tmp
    java -Djava.io.tmpdir="\$PWD/tmp" -jar /opt/raw/picard/picard.jar FixMateInformation \
      I=presplit.bam O=fixed.bam MAX_RECORDS_IN_RAM=1000000 VALIDATION_STRINGENCY=SILENT \
      AS=false SO=coordinate TMP_DIR=tmp > fixmate.log 2>&1
    java -Djava.io.tmpdir="\$PWD/tmp" -jar /opt/raw/picard/picard.jar MarkDuplicates \
      I=fixed.bam O=marked.bam M='${meta.sample_name}.metrics.txt' MAX_RECORDS_IN_RAM=1000000 \
      VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false TMP_DIR=tmp > markduplicates.log 2>&1
    samtools sort --threads ${task.cpus} -n -O bam -o '${meta.sample_name}.name_sorted.bam' marked.bam
    samtools quickcheck '${meta.sample_name}.name_sorted.bam'
    """
}
