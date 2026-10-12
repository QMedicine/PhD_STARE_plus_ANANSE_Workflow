include { ESTIMATED_COOLER_TO_H5; ESTIMATED_HIC_DEPTH_NORMALIZE; ESTIMATED_ICE_CORRECT } from '../modules/hic_estimated'

workflow HIC_ESTIMATED {
    take:
    manifest
    main:
    if (!params.hic_coolers || !params.ice_thresholds) { error 'cooler_estimated_ice requires --hic_coolers and explicit --ice_thresholds' }
    cooler_rows = channel.fromPath(params.hic_coolers,checkIfExists:true).splitCsv(header:true,sep:'\t').collect(flat:false).map { rows -> [rows] }
    cooler_inputs = manifest.combine(cooler_rows).flatMap { m,rows ->
        def states = m.samples.findAll { row -> row.assay == 'HiC' }.collect { row -> row.cell_state }.unique()
        if (rows.size()!=states.size() || rows.collect { row -> row.cell_state }.unique().sort()!=states.sort()) { error 'Supply exactly one merged cooler for each Hi-C cell_state' }
        rows.collect { row -> tuple(row.cell_state,file(row.path,checkIfExists:true)) }
    }
    thresholds = channel.fromPath(params.ice_thresholds,checkIfExists:true).splitCsv(header:true,sep:'\t').collect(flat:false).map { rows ->
        if (!rows || rows.collect { row -> row.cell_state }.unique().size()!=rows.size()) { error 'ICE threshold states must be unique' }
        rows.each { row ->
            if (row.resolution.toInteger()!=params.chosen_res || !Double.isFinite(row.lower.toDouble()) || !Double.isFinite(row.upper.toDouble()) || row.lower.toDouble()>=row.upper.toDouble() || !row.provenance) { error 'Invalid estimated ICE threshold record' }
        }
        [rows]
    }
    ESTIMATED_COOLER_TO_H5(cooler_inputs)
    group = ESTIMATED_COOLER_TO_H5.out.matrix.collect(flat:false).map { rows -> tuple(rows.collect { row -> row[0] },rows.collect { row -> row[1] }) }
    ESTIMATED_HIC_DEPTH_NORMALIZE(group)
    normalized = ESTIMATED_HIC_DEPTH_NORMALIZE.out.matrices.flatMap { states,matrices ->
        states.collect { state -> tuple(state,matrices.find { matrix -> matrix.name==state+'.normalized.h5' }) }
    }
    correction_inputs = normalized.combine(thresholds).map { state,matrix,rows ->
        def matched=rows.findAll { row -> row.cell_state==state }
        if (matched.size()!=1) { error "Missing ICE thresholds for ${state}" }
        tuple(state,matrix,matched[0].lower.toDouble(),matched[0].upper.toDouble(),matched[0].provenance)
    }
    ESTIMATED_ICE_CORRECT(correction_inputs)
    emit:
    matrix = ESTIMATED_ICE_CORRECT.out.matrix
}
