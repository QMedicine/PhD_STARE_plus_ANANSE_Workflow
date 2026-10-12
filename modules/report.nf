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
    output:
    path '*.html', emit: html
    path '*.json', emit: checksums
    script:
    """
    build_report.py --source '${source}' --target '${target}' --mode '${mode}' --ranking "${ranking}" \
      --samplesheet "${samplesheet}" --pdfs "${pdfs}" --versions "${versions}" --output '${source}_to_${target}_${mode}.html'
    """
}
