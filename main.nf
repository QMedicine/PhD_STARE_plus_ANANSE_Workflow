nextflow.enable.dsl=2
include { PHD_STARE_PLUS_ANANSE_WORKFLOW } from './workflows/phd_stare_plus_ananse_workflow'
workflow {
    if (!params.samplesheet || !params.source2target) { error '--samplesheet and --source2target are required' }
    if (!(params.stage in ['preflight','rna_recount','historical_replay','network_replay','bam_historical_hic_checkpoint','atac_historical_hic_checkpoint'])) { error "Unsupported stage: ${params.stage}" }
    if (!(params.reverse_contrast instanceof Boolean)) {
        def flag = params.reverse_contrast.toString().toLowerCase()
        if (!(flag in ['true','false'])) { error '--reverse_contrast must be TRUE or FALSE' }
        params.reverse_contrast = flag == 'true'
    }
    if (params.atac_norm != 'loess_offsets' || params.hicexp_norm != 'ICE' || params.hicexp_depth_norm != 'smallest') { error 'Alternative normalization methods are not implemented in this exact-method draft' }
    if (params.empty_plot_policy != 'fail') { error 'Only strict plotting failure is enabled pending a user-approved filter revision' }
    if (!params.atac_modes || params.atac_modes.unique().size() != params.atac_modes.size() || params.atac_modes.any { !(it in ['historical_column_order','corrected_column_order']) }) {
        error '--atac_modes must contain unique, supported ATAC modes'
    }
    if ('historical_column_order' in params.atac_modes && (!params.historical_header_samples || !params.historical_data_samples)) {
        error 'Historical ATAC mode requires explicit --historical_header_samples and --historical_data_samples'
    }
    PHD_STARE_PLUS_ANANSE_WORKFLOW(file(params.samplesheet,checkIfExists:true),file(params.source2target,checkIfExists:true))
}
