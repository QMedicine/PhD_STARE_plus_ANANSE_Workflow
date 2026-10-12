include { RNA_RECONSTRUCTION } from './rna_reconstruction'
include { VALIDATE_INPUTS } from '../modules/preflight'
include { RNA_FEATURECOUNTS; RNA_MERGE_RECOUNT; RNA_TPM_FILTER; RNA_DESEQ2 } from '../modules/rna'
include { STARE_TO_NETWORK; ANANSE_INFLUENCE } from '../modules/network'
include { PLOT_TOP20; SELF_CONTAINED_REPORT } from '../modules/report'
include { ATAC_BLACKLIST; ATAC_UNION_SUMMITS; ATAC_WINDOW_COUNTS; ATAC_ASSEMBLE; ATAC_NORMALIZE } from '../modules/atac'
include { HIC_CONTACTS_FROM_H5; STARE_GABC } from '../modules/stare'

workflow PHD_STARE_PLUS_ANANSE_WORKFLOW {
    take:
    samplesheet
    source2target

    main:
    VALIDATE_INPUTS(samplesheet, source2target)
    if (params.stage == 'rna_recount') {
        rna_manifest = VALIDATE_INPUTS.out.manifest.map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
        RNA_RECONSTRUCTION(rna_manifest,samplesheet,source2target)
    }
    if (params.stage in ['historical_replay','network_replay','bam_historical_hic_checkpoint','atac_historical_hic_checkpoint']) {
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
            if (!params.hic_checkpoints || !params.selected_counts || !params.gene_lengths) {
                error 'bam_historical_hic_checkpoint requires --hic_checkpoints, --selected_counts, --gene_lengths'
            }
            atac_samples = manifest.flatMap { m -> m.samples.findAll { it.assay == 'ATAC' }.collect { r ->
                tuple(r, file(r.bam, checkIfExists: true), file(r.peak, checkIfExists: true))
            } }
            rna_samples = manifest.flatMap { m -> m.samples.findAll { it.assay == 'RNA' }.collect { r ->
                tuple(r, file(r.bam, checkIfExists: true))
            } }
            if (params.stage == 'atac_historical_hic_checkpoint') {
                tpm = channel.value(file(params.tpm,checkIfExists:true))
                de_ch = channel.fromPath(params.de_manifest,checkIfExists:true).splitCsv(header:true,sep:'\t').map { row -> tuple(row.source,row.target,file(row.path,checkIfExists:true)) }
            } else {
                RNA_FEATURECOUNTS(rna_samples, file(params.gtf, checkIfExists: true))
                RNA_MERGE_RECOUNT(RNA_FEATURECOUNTS.out.counts.map { meta, counts -> counts }.collect(), file(params.selected_counts, checkIfExists: true))
                RNA_TPM_FILTER(RNA_MERGE_RECOUNT.out.counts, file(params.gene_lengths, checkIfExists: true), file(params.idmap, checkIfExists: true))
                tpm = RNA_TPM_FILTER.out.tpm.first()
                forward_transitions = Channel.fromPath(params.source2target, checkIfExists: true).splitCsv(header: true, sep: '\t').map { tuple(it.source, it.target) }
                RNA_DESEQ2(forward_transitions, RNA_TPM_FILTER.out.counts.first(), samplesheet)
                de_ch = params.reverse_contrast ? RNA_DESEQ2.out.forward.mix(RNA_DESEQ2.out.reverse) : RNA_DESEQ2.out.forward
            }
            ATAC_BLACKLIST(atac_samples, file(params.atac_blacklist, checkIfExists: true))
            ATAC_UNION_SUMMITS(ATAC_BLACKLIST.out.bam.map { meta, bam, summit -> summit }.collect(), file(params.chrom_sizes, checkIfExists: true), file(params.summit_blacklist, checkIfExists: true))
            ATAC_WINDOW_COUNTS(ATAC_BLACKLIST.out.bam, ATAC_UNION_SUMMITS.out.windows.first())
            ATAC_ASSEMBLE(Channel.fromList(params.atac_modes), ATAC_WINDOW_COUNTS.out.counts.collect(), samplesheet)
            ATAC_NORMALIZE(ATAC_ASSEMBLE.out.counts)
            hic_matrices = Channel.fromPath(params.hic_checkpoints, checkIfExists: true).splitCsv(header: true, sep: '\t').map { row -> tuple(row.cell_state, file(row.path, checkIfExists: true)) }
            HIC_CONTACTS_FROM_H5(hic_matrices)
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
            stare_inputs = state_columns.combine(ATAC_NORMALIZE.out.activity, by: 0).map { mode,state,columns,activity -> tuple(state,mode,columns,activity) }.combine(HIC_CONTACTS_FROM_H5.out.contacts, by: 0).map { state,mode,columns,activity,contacts -> tuple(mode,state,columns,activity,contacts) }
            STARE_GABC(stare_inputs, samplesheet, file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), file(params.psem, checkIfExists: true), file("${projectDir}/historical/scripts/06_stare/modified_STARE.sh"))
            affinities = STARE_GABC.out.affinities
        }
        if (params.stage != 'network_replay') {
        STARE_TO_NETWORK(affinities, tpm, file(params.idmap, checkIfExists: true), file(params.gtf, checkIfExists: true), file("${projectDir}/historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py"), file("${projectDir}/historical/scripts/07_network_conversion/02_convert_ensembl_to_symbol.py"))
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
        plot_inputs = ANANSE_INFLUENCE.out.results.map { mode,source,target,ranking,diffnet,de ->
            tuple(mode,source,target,ranking,diffnet,de,file(params.plot_template ?: "${projectDir}/bin/plot_branching_template.R",checkIfExists:true))
        }
        PLOT_TOP20(plot_inputs)
        SELF_CONTAINED_REPORT(PLOT_TOP20.out.pdfs, samplesheet)
        }
    }
}
