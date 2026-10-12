process PLOT_TOP20 {
    tag "${mode}:${source}->${target}"
    label 'r_plot'
    cpus 4
    memory '32 GB'
    time '8h'
    publishDir { "${params.outdir}/${mode}/${source}_to_${target}/plots" }, mode: 'copy'
    input:
    tuple val(mode), val(source), val(target), path(ranking), path(diffnetwork), path(de), path(notebook)
    output:
    tuple val(mode), val(source), val(target), path(ranking), path('plots'), path('versions.txt'), emit: pdfs
    script:
    def display = params.plot_settings["${source}->${target}"] ?: [limits:params.plot_fc_limits,breaks:params.plot_fc_breaks]
    """
    plot_top20.R "${notebook}" "${ranking}" "${diffnetwork}" "${de}" '${source}' '${target}' plots ${display.limits[0]} ${display.limits[1]} '${display.breaks.join(',')}' '${params.empty_plot_policy}'
    """
}

process SELF_CONTAINED_REPORT {
    tag "${mode}:${source}->${target}"
    label 'ananse'
    publishDir { "${params.outdir}/${mode}/${source}_to_${target}/reports" }, mode: 'copy'
    input:
    tuple val(mode), val(source), val(target), path(ranking), path(pdfs), path(versions)
    path samplesheet
    path software_lock
    output:
    path '*.html', emit: html
    path '*.json', emit: checksums
    script:
    def keys = ['stage','process_fastq','rna_counts_mode','chosen_res','summit_flank','stare_window','enhancer_window','gabc_threshold','min_binding','influence_edges','min_count','min_count_tf','min_tpm','min_tpm_tf','union_windows','reverse_contrast','gtf','fasta','psem','atac_blacklist','hic_checkpoints']
    def context = keys.collectEntries { key -> [(key):params[key]] }
    context.containers = ['ananse','r_plot','tpm','featurecounts','hicexplorer','atac_tools','bedtools','stare','deseq2','raw_tools'].collectEntries { role -> [(role):params[role+'_container']] }
    def encoded = groovy.json.JsonOutput.toJson(context).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > methods_context.json
    build_report.py --source '${source}' --target '${target}' --mode '${mode}' --ranking "${ranking}" \
      --samplesheet "${samplesheet}" --pdfs "${pdfs}" --versions "${versions}" --context methods_context.json \
      --software-lock "${software_lock}" --output '${source}_to_${target}_${mode}.html'
    """
}
