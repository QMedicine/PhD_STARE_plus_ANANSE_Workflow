nextflow.enable.dsl=2
include { PREPARE_TF_EXPRESSION; PLOT_TOP20; SELF_CONTAINED_REPORT; report_variants } from './modules/report'
// Report-only entry point consumes explicit, previously produced score tables.
workflow {
    if (!params.report_manifest || !params.samplesheet) { error '--report_manifest and --samplesheet are required' }
    def samples = file(params.samplesheet,checkIfExists:true)
    def lines = samples.readLines()
    def stateIndex = lines[0].split('\t',-1).toList().indexOf('cell_state')
    if (stateIndex < 0) { error 'samplesheet requires cell_state' }
    def states = lines.drop(1).findAll { line -> line.trim() }.collect { line -> line.split('\t',-1)[stateIndex] }.unique()
    inputs = channel.fromPath(params.report_manifest,checkIfExists:true).splitCsv(header:true,sep:'\t').flatMap { row ->
        if ([row.mode,row.source,row.target].any { value -> !(value ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) }) { error 'Report mode and states require safe label tokens' }
        if ([row.ranking,row.diffnetwork,row.de].any { value -> !value.startsWith('/') || value.find(/[\n\r\t'"`\$]/) }) { error 'Report input paths must be safe absolute paths' }
        if (!(row.source in states) || !(row.target in states) || row.source==row.target) { error 'Report states must match distinct samplesheet cell_state values' }
        report_variants(row.mode,row.source,row.target,file(row.ranking,checkIfExists:true),file(row.diffnetwork,checkIfExists:true),file(row.de,checkIfExists:true),
            file("${projectDir}/bin/plot_branching_template.R",checkIfExists:true),file("${projectDir}/bin/plot_top10_template.R",checkIfExists:true))
    }
    report_font = params.heatmap_arial_font ? [file(params.heatmap_arial_font,checkIfExists:true)] : []
    PREPARE_TF_EXPRESSION(inputs, file(params.tpm,checkIfExists:true), samples, report_font)
    PLOT_TOP20(PREPARE_TF_EXPRESSION.out.prepared, report_font)
    SELF_CONTAINED_REPORT(PLOT_TOP20.out.pdfs,samples,file("${projectDir}/containers/software_lock_v01.json",checkIfExists:true))
}
