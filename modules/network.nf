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
    path idmap
    path gtf
    path converter
    path symbol_converter
    output:
    tuple val(mode), val(state), path("${state}_SYMBOLS.tsv"), emit: network
    path 'raw/network_meta.json', emit: metadata
    script:
    """
    python "${converter}" --stare ${affinities.collect { '"'+it+'"' }.join(' ')} --tpm "${tpm}" \
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
    tuple val(mode), val(source), val(target), path('influence.txt'), path('influence_diffnetwork.tsv'), path(de), emit: results
    script:
    """
    export XDG_CACHE_HOME="\$PWD/cache"
    mkdir -p "\$XDG_CACHE_HOME"
    ananse influence -s "${source_network}" -t "${target_network}" -d "${de}" -a "${gtf}" \
      -i ${params.influence_edges} -o influence.txt -n ${task.cpus}
    """
}
