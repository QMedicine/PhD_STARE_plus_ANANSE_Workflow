// Explicit exploratory route: retained pre-blacklist coolers and estimated ICE thresholds.
process ESTIMATED_COOLER_TO_H5 {
    tag state
    label 'hicexplorer'
    cpus 2
    memory '24 GB'
    time '6h'
    input:
    tuple val(state), path(cooler)
    output:
    tuple val(state), path("${state}.uncorrected.h5"), emit: matrix
    script:
    """
    hicConvertFormat --matrices "${cooler}" --inputFormat cool --outputFormat h5 --outFileName '${state}.uncorrected.h5'
    validate_hic_resolution.py --matrix '${state}.uncorrected.h5' --resolution ${params.chosen_res}
    """
}

process ESTIMATED_HIC_DEPTH_NORMALIZE {
    label 'hicexplorer'
    cpus 2
    memory '48 GB'
    time '6h'
    input:
    tuple val(states), path(matrices)
    output:
    tuple val(states), path('normalized/*.h5'), emit: matrices
    script:
    def inputs = matrices.collect { matrix -> '"'+matrix+'"' }.join(' ')
    def outputs = states.collect { state -> "'normalized/${state}.normalized.h5'" }.join(' ')
    """
    mkdir normalized
    hicNormalize --matrices ${inputs} --normalize ${params.hicexp_depth_norm} --outFileName ${outputs}
    """
}

process ESTIMATED_ICE_CORRECT {
    tag state
    label 'hicexplorer'
    cpus 2
    memory '48 GB'
    time '12h'
    publishDir { "${params.outdir}/estimated_ice/${state}" }, mode: 'copy'
    input:
    tuple val(state), path(matrix), val(lower), val(upper), val(provenance)
    output:
    tuple val(state), path("${state}.estimated_ICE.h5"), emit: matrix
    path '*.diagnostic.pdf', emit: diagnostic
    path 'hic_estimated_provenance.json', emit: provenance
    path 'versions.txt', emit: versions
    script:
    def receipt = groovy.json.JsonOutput.toJson([cell_state:state,chosen_res:params.chosen_res,lower:lower,upper:upper,threshold_provenance:provenance,status:'ESTIMATED_PRE_BLACKLIST_NOT_HISTORICALLY_AUTHENTICATED',depth_normalization:params.hicexp_depth_norm,correction:params.hicexp_norm]).bytes.encodeBase64().toString()
    """
    hicCorrectMatrix diagnostic_plot --matrix "${matrix}" --plotName '${state}.diagnostic.pdf'
    hicCorrectMatrix correct --matrix "${matrix}" --correctionMethod ${params.hicexp_norm} \
      --filterThreshold ${lower} ${upper} --outFileName '${state}.estimated_ICE.h5'
    validate_hic_resolution.py --matrix '${state}.estimated_ICE.h5' --resolution ${params.chosen_res}
    printf '%s' '${receipt}' | base64 -d > hic_estimated_provenance.json
    hicCorrectMatrix --version > versions.txt 2>&1
    """
}
