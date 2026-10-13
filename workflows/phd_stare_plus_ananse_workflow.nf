include { RAW_PREPARE } from './raw_prepare'
include { HIC_ESTIMATED } from './hic_estimated'
include { DIFFHIC_COUNTS } from './diffhic_counts'
include { RNA_RECONSTRUCTION } from './rna_reconstruction'
include { VALIDATE_INPUTS } from '../modules/preflight'
include { RNA_FEATURECOUNTS; RNA_MERGE_RECOUNT; RNA_ANNOTATION_COUNTS; RNA_TPM_FILTER; RNA_DESEQ2 } from '../modules/rna'
include { STARE_TO_NETWORK; ANANSE_INFLUENCE } from '../modules/network'
include { PLOT_TOP20; SELF_CONTAINED_REPORT; report_variants } from '../modules/report'
include { ATAC_BLACKLIST; ATAC_UNION_SUMMITS; ATAC_WINDOW_COUNTS; ATAC_ASSEMBLE; ATAC_NORMALIZE } from '../modules/atac'
include { HIC_CONTACTS_FROM_H5; STARE_GABC } from '../modules/stare'

workflow PHD_STARE_PLUS_ANANSE_WORKFLOW {
    take:
    samplesheet
    source2target

    main:
    VALIDATE_INPUTS(samplesheet, source2target)
    if (params.stage == 'hic_diffhic_counts') {
        diffhic_manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        if (params.process_fastq) { RAW_PREPARE(diffhic_manifest) }
        DIFFHIC_COUNTS(diffhic_manifest,samplesheet,params.process_fastq ? RAW_PREPARE.out.hic_bam : channel.empty())
    }
    if (params.stage == 'cooler_estimated_ice') {
        hic_manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        HIC_ESTIMATED(hic_manifest)
    }
    if (params.stage == 'raw_prepare') {
        raw_manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        RAW_PREPARE(raw_manifest)
    }
    if (params.stage == 'rna_recount') {
        rna_manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        if (params.process_fastq) { RAW_PREPARE(rna_manifest) }
        RNA_RECONSTRUCTION(rna_manifest,samplesheet,source2target,params.process_fastq ? RAW_PREPARE.out.rna_bam : channel.empty())
    }
    if (params.stage in ['historical_replay','network_replay','bam_historical_hic_checkpoint','atac_historical_hic_checkpoint','atac_diffhic_hic_contacts','bam_diffhic_hic_contacts']) {
        manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        transitions = manifest.flatMap { it.transitions.collect { t -> tuple(t.source, t.target) } }
        if (params.stage == 'network_replay') {
            if (!params.network_manifest || !params.de_manifest) { error 'network_replay requires --network_manifest and --de_manifest' }
            replay_networks = channel.fromPath(params.network_manifest,checkIfExists:true).splitCsv(header:true,sep:'\t').map { row -> tuple(row.mode,row.cell_state,file(row.path,checkIfExists:true)) }
            de_ch = channel.fromPath(params.de_manifest,checkIfExists:true).splitCsv(header:true,sep:'\t').map { row -> tuple(row.source,row.target,file(row.path,checkIfExists:true)) }
        } else if (params.stage == 'historical_replay') {
            if (!params.affinities || !params.de_manifest || !params.tpm) {
                error 'historical_replay requires --affinities, --de_manifest, --tpm'
            }
            affinities = Channel.fromPath(params.affinities, checkIfExists: true).splitCsv(header: true, sep: '\t').map { row ->
                tuple(row.mode ?: 'archived_affinity_replay', row.cell_state, file(row.path, checkIfExists: true))
            }.groupTuple(by: [0,1])
            tpm = Channel.value(file(params.tpm, checkIfExists: true))
            de_ch = Channel.fromPath(params.de_manifest, checkIfExists: true).splitCsv(header: true, sep: '\t').map { row ->
                tuple(row.source, row.target, file(row.path, checkIfExists: true))
            }
        } else {
            if (!params.hic_checkpoints && !params.hic_contacts_from_diffhic) {
                error 'Hi-C checkpoint routes require --hic_checkpoints'
            }
            if (params.process_fastq) { RAW_PREPARE(manifest) }
            atac_samples = manifest.flatMap { m -> m.samples.findAll { it.assay == 'ATAC' && it.bam }.collect { r ->
                tuple(r, file(r.bam, checkIfExists: true), file(r.peak, checkIfExists: true))
            } }
            rna_samples = manifest.flatMap { m -> m.samples.findAll { it.assay == 'RNA' && it.bam }.collect { r ->
                tuple(r, file(r.bam, checkIfExists: true))
            } }
            if (params.process_fastq) {
                rna_samples = rna_samples.mix(RAW_PREPARE.out.rna_bam)
                // Raw ATAC output has already undergone the historical blacklist step.
                raw_atac = RAW_PREPARE.out.atac_bam_peaks
            }
            if (params.stage in ['atac_historical_hic_checkpoint','atac_diffhic_hic_contacts']) {
                tpm = channel.value(file(params.tpm,checkIfExists:true))
                de_ch = channel.fromPath(params.de_manifest,checkIfExists:true).splitCsv(header:true,sep:'\t').map { row -> tuple(row.source,row.target,file(row.path,checkIfExists:true)) }
            } else {
                RNA_FEATURECOUNTS(rna_samples, file(params.gtf, checkIfExists: true))
                if (params.rna_counts_mode == 'historical_selection') {
                    RNA_MERGE_RECOUNT(RNA_FEATURECOUNTS.out.counts.map { _meta, counts -> counts }.collect(), file(params.selected_counts, checkIfExists: true))
                    expression_counts = RNA_MERGE_RECOUNT.out.counts
                    expression_lengths = channel.value(file(params.gene_lengths,checkIfExists:true))
                    expression_map = channel.value(file(params.idmap,checkIfExists:true))
                } else {
                    RNA_ANNOTATION_COUNTS(RNA_FEATURECOUNTS.out.counts.map { _meta, counts -> counts }.collect(),samplesheet)
                    expression_counts = RNA_ANNOTATION_COUNTS.out.counts
                    expression_lengths = params.gene_lengths ? channel.value(file(params.gene_lengths,checkIfExists:true)) : RNA_ANNOTATION_COUNTS.out.lengths.first()
                    expression_map = params.idmap ? channel.value(file(params.idmap,checkIfExists:true)) : RNA_ANNOTATION_COUNTS.out.idmap.first()
                }
                RNA_TPM_FILTER(expression_counts, expression_lengths, expression_map)
                tpm = RNA_TPM_FILTER.out.tpm.first()
                forward_transitions = manifest.flatMap { m -> m.primary_transitions.collect { r -> tuple(r.source,r.target) } }
                RNA_DESEQ2(forward_transitions, RNA_TPM_FILTER.out.counts.first(), samplesheet)
                de_ch = params.reverse_contrast ? RNA_DESEQ2.out.forward.mix(RNA_DESEQ2.out.reverse) : RNA_DESEQ2.out.forward
            }
            ATAC_BLACKLIST(atac_samples, file(params.atac_blacklist, checkIfExists: true))
            prepared_atac = params.process_fastq ? ATAC_BLACKLIST.out.bam.mix(raw_atac) : ATAC_BLACKLIST.out.bam
            if (params.union_windows) {
                union_windows_ch = channel.value(file(params.union_windows,checkIfExists:true))
            } else {
                ATAC_UNION_SUMMITS(prepared_atac.map { _meta,_bam,summit -> summit }.collect(),file(params.chrom_sizes,checkIfExists:true),file(params.summit_blacklist,checkIfExists:true))
                union_windows_ch = ATAC_UNION_SUMMITS.out.windows.first()
            }
            ATAC_WINDOW_COUNTS(prepared_atac,union_windows_ch)
            ATAC_ASSEMBLE(Channel.fromList(params.atac_modes), ATAC_WINDOW_COUNTS.out.counts.collect(), samplesheet)
            ATAC_NORMALIZE(ATAC_ASSEMBLE.out.counts)
            if (params.hic_contacts_from_diffhic) {
                DIFFHIC_COUNTS(manifest,samplesheet,params.process_fastq ? RAW_PREPARE.out.hic_bam : channel.empty())
                hic_contacts_ch = DIFFHIC_COUNTS.out.contacts
            } else {
                hic_matrices = Channel.fromPath(params.hic_checkpoints, checkIfExists: true).splitCsv(header: true, sep: '\t').map { row -> tuple(row.cell_state, file(row.path, checkIfExists: true)) }
                HIC_CONTACTS_FROM_H5(hic_matrices)
                hic_contacts_ch = HIC_CONTACTS_FROM_H5.out.contacts
            }
            state_columns = manifest.flatMap { m ->
                def states = m.samples.findAll { it.assay == 'ATAC' }.collect { it.cell_state }.unique()
                states.collectMany { state ->
                    params.atac_modes.collect { mode ->
                        def rows = m.samples.findAll { it.assay == 'ATAC' }
                        if (mode == 'historical_column_order') {
                            rows = params.historical_header_samples.collect { sample ->
                                def matching = rows.findAll { it.sample_name == sample }
                                if (matching.size() != 1) { error 'Historical header samples must match ATAC identities exactly' }
                                matching[0]
                            }
                        }
                        def columns = rows.withIndex().findAll { r,i -> r.cell_state == state }.collect { r,i -> i+5 }.join(',')
                        tuple(mode,state,columns)
                    }
                }
            }
            stare_inputs = state_columns.combine(ATAC_NORMALIZE.out.activity, by: 0).map { mode,state,columns,activity -> tuple(state,mode,columns,activity) }.combine(hic_contacts_ch, by: 0).map { state,mode,columns,activity,contacts -> tuple(params.hic_contacts_from_diffhic ? 'diffhic_loess_fixed_grid__'+mode : mode,state,columns,activity,contacts) }
            STARE_GABC(stare_inputs, samplesheet, file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), file(params.psem, checkIfExists: true), file("${projectDir}/historical/scripts/06_stare/modified_STARE.sh"))
            affinities = STARE_GABC.out.affinities
        }
        if (params.stage != 'network_replay') {
        conversion_map = params.stage in ['bam_historical_hic_checkpoint','bam_diffhic_hic_contacts'] ? expression_map : channel.value(file(params.idmap,checkIfExists:true))
        STARE_TO_NETWORK(affinities, tpm, conversion_map, file(params.gtf, checkIfExists: true), file("${projectDir}/historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py"), file("${projectDir}/historical/scripts/07_network_conversion/02_convert_ensembl_to_symbol.py"))
        // Value-channel collection is bounded by the number of states; explicit keys drive lookup.
        }
        network_ch = params.stage == 'network_replay' ? replay_networks : STARE_TO_NETWORK.out.network
        networks = network_ch.collect(flat:false).map { nets -> [nets] }
        influence_inputs = transitions.combine(networks).combine(de_ch.collect(flat:false).map { des -> [des] }).flatMap { source,target,nets,des ->
            def selected_de = des.findAll { it[0] == source && it[1] == target }
            if (selected_de.size() != 1) { error "Expected one DE file for ${source}->${target}" }
            def modes = nets.collect { it[0] }.unique()
            modes.collect { mode ->
                def src = nets.findAll { it[0] == mode && it[1] == source }
                def tgt = nets.findAll { it[0] == mode && it[1] == target }
                if (src.size() != 1 || tgt.size() != 1) { error "Expected source and target networks for ${mode}:${source}->${target}" }
                tuple(mode,source,target,src[0][2],tgt[0][2],selected_de[0][2])
            }
        }
        ANANSE_INFLUENCE(influence_inputs, file(params.gtf, checkIfExists: true))
        if (params.make_reports) {
        plot_inputs = ANANSE_INFLUENCE.out.results.flatMap { mode,source,target,ranking,diffnet,de ->
            report_variants(mode,source,target,ranking,diffnet,de,file(params.plot_template ?: "${projectDir}/bin/plot_branching_template.R",checkIfExists:true),file("${projectDir}/bin/plot_top10_template.R",checkIfExists:true))
        }
        PLOT_TOP20(plot_inputs)
        SELF_CONTAINED_REPORT(PLOT_TOP20.out.pdfs, samplesheet, file("${projectDir}/containers/software_lock_v01.json",checkIfExists:true))
        }
    }
}
