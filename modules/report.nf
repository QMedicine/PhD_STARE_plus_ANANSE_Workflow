// Explicit display variants share the same ranking and directional DE bytes.
def report_variants(mode, source, target, ranking, diffnet, de, notebook, top10) {
    def settings = [limits:params.plot_fc_limits, breaks:params.plot_fc_breaks,
        min_edge_weight_plot:params.plot_min_edge_weight, padj_cutoff:params.plot_padj_cutoff,
        abs_lfc_cutoff:params.plot_abs_lfc_cutoff, empty_plot_policy:params.empty_plot_policy]
    settings.putAll(params.plot_settings["${source}->${target}"] ?: [:])
    def variants = [tuple(mode, source, target, 'main', settings, ranking, diffnet, de, notebook, top10)]
    if (params.make_relaxed_reports) {
        def relaxed = settings + [min_edge_weight_plot:settings.min_edge_weight_plot * params.relaxed_plot_edge_factor,
            abs_lfc_cutoff:Math.min(settings.abs_lfc_cutoff as double, params.relaxed_plot_abs_lfc_cutoff as double)]
        def label = "alternate_filters_edge_${relaxed.min_edge_weight_plot}_padj_${relaxed.padj_cutoff}_absLFC_${relaxed.abs_lfc_cutoff}"
        variants.add(tuple(mode, source, target, label, relaxed, ranking, diffnet, de, notebook, top10))
    }
    variants
}

process PREPARE_TF_EXPRESSION {
    tag "${mode}:${source}->${target}:${variant}"
    label 'ananse'
    cpus 1
    memory '4 GB'
    time '15m'
    input:
    tuple val(mode), val(source), val(target), val(variant), val(settings), path(ranking), path(diffnetwork), path(de, stageAs:'directional_de.tsv'), path(notebook), path(top10)
    path tpm
    path samplesheet
    path arial_font
    output:
    tuple val(mode), val(source), val(target), val(variant), val(settings), path(ranking), path(diffnetwork), path(de), path(notebook), path(top10), path('top20_tf_expression.tsv'), path('heatmap_settings.json'), emit: prepared
    script:
    def fontPath = arial_font instanceof List ? (arial_font ? arial_font[0] : null) : arial_font
    def fontOption = fontPath ? '--font-file "'+fontPath+'"' : ''
    """
    prepare_tf_expression.py --ranking "${ranking}" --de "${de}" --tpm "${tpm}" \
      --samplesheet "${samplesheet}" --source '${source}' --target '${target}' \
      --output top20_tf_expression.tsv --settings heatmap_settings.json \
      --tpm-limits ${params.heatmap_tpm_limits.join(' ')} --fc-limits ${params.heatmap_fc_limits.join(' ')} ${fontOption}
    """
}

process PLOT_TOP20 {
    tag "${mode}:${source}->${target}:${variant}"
    label 'r_plot'
    cpus 4
    memory '32 GB'
    time '8h'
    publishDir { "${params.outdir}/${mode}/${source}_to_${target}/${variant}" }, mode: 'copy'
    input:
    tuple val(mode), val(source), val(target), val(variant), val(settings), path(ranking), path(diffnetwork), path(de, stageAs:'directional_de.tsv'), path(notebook), path(top10), path(expression_table), path(heatmap_settings)
    path arial_font
    output:
    tuple val(mode), val(source), val(target), val(variant), path(ranking), path('plots'), path('versions.txt'), emit: pdfs
    script:
    def encoded = groovy.json.JsonOutput.toJson(settings).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > plot_settings.json
    mkdir -p plots
    plot_tf_expression.R "${expression_table}" "${heatmap_settings}" plots/TF_rank_RNA_expression.pdf
    cp "${expression_table}" "${heatmap_settings}" plots/
    plot_top20.R "${notebook}" "${ranking}" "${diffnetwork}" "${de}" '${source}' '${target}' plots plot_settings.json "${top10}"
    """
}

process SELF_CONTAINED_REPORT {
    tag "${mode}:${source}->${target}:${variant}"
    label 'ananse'
    publishDir { "${params.outdir}/${mode}/${source}_to_${target}/${variant}/reports" }, mode: 'copy'
    input:
    tuple val(mode), val(source), val(target), val(variant), path(ranking), path(pdfs), path(versions)
    path samplesheet
    path software_lock
    output:
    path '*.html', emit: html
    path '*.json', emit: checksums
    script:
    def keys = ['stage','process_fastq','rna_counts_mode','chosen_res','summit_flank','stare_window','enhancer_window','gabc_threshold','min_binding','influence_edges','min_count','min_count_tf','min_tpm','min_tpm_tf','union_windows','reverse_contrast','gtf','fasta','psem','atac_blacklist','hic_checkpoints','hic_contacts_from_diffhic','hic_mapq','hic_chim_dist','hic_max_frag','hic_min_inward','hic_min_outward','hic_restriction_sites','hic_overhang','diffhic_norm','report_manifest','hic_alignment_mode','hic_ligation_signature','tpm_expression_mode','heatmap_tpm_limits','heatmap_fc_limits','heatmap_arial_font']
    def context = keys.collectEntries { key -> [(key):params[key]] }
    context.plot_variant = variant
    context.containers = ['ananse','r_plot','tpm','featurecounts','hicexplorer','atac_tools','bedtools','stare','deseq2','raw_tools','hic_raw','diffhic'].collectEntries { role -> [(role):params[role+'_container']] }
    def encoded = groovy.json.JsonOutput.toJson(context).bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded}' | base64 -d > methods_context.json
    build_report.py --source '${source}' --target '${target}' --mode '${mode}' --ranking "${ranking}" \
      --samplesheet "${samplesheet}" --pdfs "${pdfs}" --versions "${versions}" --context methods_context.json \
      --software-lock "${software_lock}" --output '${source}_to_${target}_${mode}_${variant}.html'
    """
}
