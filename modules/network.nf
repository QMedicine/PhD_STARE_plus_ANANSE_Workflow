process STARE_TO_NETWORK {
    tag "${mode}:${state}"
    label 'ananse'
    cpus 4
    memory '64 GB'
    time '24h'
    publishDir { "${params.outdir}/${mode}/networks/${state}" }, mode: 'copy'
    input:
    tuple val(mode), val(state), path(affinities, stageAs: 'affinities/*')
    path tpm
    path samplesheet
    path idmap
    path gtf
    path converter
    path symbol_converter
    output:
    tuple val(mode), val(state), path("${state}_SYMBOLS.tsv"), emit: network
    path 'raw/network_meta.json', emit: metadata
    tuple val(mode), val(state), path('tpm_selection.json'), emit: expression_selection
    script:
    """
    select_state_tpm.py --tpm "${tpm}" --samplesheet "${samplesheet}" --state '${state}' \
      --mode '${params.tpm_expression_mode}' --output selected_state.tpm.tsv --metadata tpm_selection.json
    python "${converter}" --stare ${affinities.collect { '"'+it+'"' }.join(' ')} --tpm selected_state.tpm.tsv \
      --idmap "${idmap}" --outdir raw --basename '${state}' --min-binding ${params.min_binding}
    python "${symbol_converter}" --gtf "${gtf}" --out-suffix _SYMBOLS 'raw/${state}.tsv'
    cp 'raw/${state}_SYMBOLS.tsv' '${state}_SYMBOLS.tsv'
    """
}

process ANANSE_INFLUENCE {
    tag "${mode}:${source}->${target}"
    label 'ananse'
    cpus 16
    memory '64 GB'
    time '4h'
    publishDir { "${params.outdir}/${mode}/${source}_to_${target}/tables" }, mode: 'copy'
    input:
    tuple val(mode), val(source), val(target), path(source_network), path(target_network), path(de, stageAs: 'directed_de.tsv')
    path gtf
    output:
    path 'de_compatibility.json', emit: de_compatibility
    tuple val(mode), val(source), val(target), path('influence.txt'), path('influence_diffnetwork.tsv'), path(de), emit: results
    script:
    """
    export XDG_CACHE_HOME="\$PWD/cache"
    mkdir -p "\$XDG_CACHE_HOME"
    prepare_ananse_de.py --input "${de}" --output ananse_de.tsv
    ananse influence -s "${source_network}" -t "${target_network}" -d ananse_de.tsv -a "${gtf}" \
      -i ${params.influence_edges} -o influence.txt -n ${task.cpus}
    """
}
