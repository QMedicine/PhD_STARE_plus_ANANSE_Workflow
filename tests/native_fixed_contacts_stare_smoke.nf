nextflow.enable.dsl=2
include { DIFFHIC_FRAGMENTS; DIFFHIC_PREPARE_PAIRS; DIFFHIC_FIXED_CONTACTS } from '../modules/diffhic'
include { STARE_GABC } from '../modules/stare'
workflow {
    manifest = channel.fromPath(params.fixture_manifest,checkIfExists:true).map { json -> new groovy.json.JsonSlurper().parseText(json.text) }
    samples = manifest.flatMap { m -> m.samples.findAll { r -> r.assay=='HiC' }.collect { r -> tuple(r,file(r.bam,checkIfExists:true)) } }
    DIFFHIC_FRAGMENTS(file(params.fasta,checkIfExists:true))
    DIFFHIC_PREPARE_PAIRS(samples,DIFFHIC_FRAGMENTS.out.fragments.first())
    DIFFHIC_FIXED_CONTACTS(DIFFHIC_PREPARE_PAIRS.out.pairs.map { _meta,pairs -> pairs }.collect(),DIFFHIC_FRAGMENTS.out.fragments.first(),file(params.samplesheet,checkIfExists:true),file("${projectDir}/../bin/fixed_grid_helpers.R"))
    state_contacts = DIFFHIC_FIXED_CONTACTS.out.contacts.combine(manifest).flatMap { folder,m -> m.samples.findAll { r -> r.assay=='HiC' }.collect { r -> r.cell_state }.unique().collect { state -> tuple(state,folder.resolve(state)) } }
    inputs = state_contacts.map { state,contacts -> tuple('diffhic_loess_fixed_grid__synthetic_smoke',state,'5,6',file(params.fixture_activity,checkIfExists:true),contacts) }
    STARE_GABC(inputs,file(params.samplesheet,checkIfExists:true),file(params.fasta,checkIfExists:true),file(params.gtf,checkIfExists:true),file(params.psem,checkIfExists:true),file("${projectDir}/../historical/scripts/06_stare/modified_STARE.sh"))
}
