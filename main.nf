nextflow.enable.dsl=2
include { PHD_STARE_PLUS_ANANSE_WORKFLOW } from './workflows/phd_stare_plus_ananse_workflow'
workflow {
    if (!params.samplesheet || !params.source2target) { error '--samplesheet and --source2target are required' }
    if (!(params.tpm_expression_mode in ['state_specific','historical_pooled'])) { error 'Unsupported --tpm_expression_mode' }
    if (!(params.stage in ['preflight','raw_prepare','hic_diffhic_counts','atac_diffhic_hic_contacts','bam_diffhic_hic_contacts','cooler_estimated_ice','rna_recount','historical_replay','network_replay','bam_historical_hic_checkpoint','atac_historical_hic_checkpoint'])) { error "Unsupported stage: ${params.stage}" }
    ['reverse_contrast','process_fastq','make_reports','make_relaxed_reports','hic_contacts_from_diffhic'].each { name ->
        if (!(params[name] instanceof Boolean)) {
            def flag = params[name].toString().toLowerCase()
            if (!(flag in ['true','false'])) { error "--${name} must be TRUE or FALSE" }
            params[name] = flag == 'true'
        }
    }
    if (params.stage in ['atac_diffhic_hic_contacts','bam_diffhic_hic_contacts']) { params.hic_contacts_from_diffhic = true }
    if (params.atac_norm != 'loess_offsets' || params.diffhic_norm != 'loess_offsets' || params.hicexp_norm != 'ICE' || params.hicexp_depth_norm != 'smallest') { error 'Alternative normalization methods are not implemented in this exact-method draft' }
    if (params.process_fastq && params.raw_sample_names && params.stage != 'raw_prepare') { error 'Partial FASTQ sample selection is permitted only in the raw_prepare pilot stage' }
    if (params.empty_plot_policy != 'center_warning') { error 'empty_plot_policy must be center_warning for the selected plotting design' }
    if (params.hic_alignment_mode && params.hic_alignment_mode != 'diffhic_presplit_standard_compat') { error 'Unsupported Hi-C alignment mode; no implicit mapper substitution' }
    if (params.hic_alignment_mode && (params.hic_restriction_sites.size()!=1 || !(params.hic_restriction_sites[0] ==~ /[ACGT]+/) || !(params.hic_ligation_signature ==~ /[ACGT]+/) || params.hic_ligation_signature.size()%2!=0)) { error 'Standard presplit route requires one restriction site and one even-length ligation signature; cocktails require a separately qualified mapping route' }
    if (!(params.rna_counts_mode in ['historical_selection','annotation_symbol_sum'])) { error 'Unsupported --rna_counts_mode' }
    if (!(params.atac_samtools in ['samtools120','samtools']) || params.atac_sort_threads<1) { error 'Unsupported ATAC samtools selector or sort-thread count' }
    if (!(params.atac_fastq_header_mode in ['preserve','sra_spot_length'])) { error 'Unsupported --atac_fastq_header_mode' }
    if (params.stage in ['rna_recount','bam_historical_hic_checkpoint','bam_diffhic_hic_contacts'] && params.rna_counts_mode=='historical_selection' && (!params.selected_counts || !params.idmap || !params.gene_lengths)) {
        error 'historical_selection requires --selected_counts, --idmap and --gene_lengths; new datasets can explicitly select --rna_counts_mode annotation_symbol_sum'
    }
    if (!params.atac_modes || params.atac_modes.unique().size() != params.atac_modes.size() || params.atac_modes.any { !(it in ['historical_column_order','corrected_column_order']) }) {
        error '--atac_modes must contain unique, supported ATAC modes'
    }
    if ('historical_column_order' in params.atac_modes && (!params.historical_header_samples || !params.historical_data_samples)) {
        error 'Historical ATAC mode requires explicit --historical_header_samples and --historical_data_samples'
    }
    PHD_STARE_PLUS_ANANSE_WORKFLOW(file(params.samplesheet,checkIfExists:true),file(params.source2target,checkIfExists:true))
}
