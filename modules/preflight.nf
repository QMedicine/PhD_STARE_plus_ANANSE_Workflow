process VALIDATE_INPUTS {
    label 'ananse'
    publishDir "${params.outdir}/provenance", mode: 'copy'
    input:
    path samplesheet
    path source2target
    output:
    path 'validated_inputs.json', emit: manifest
    script:
    def reverse = params.reverse_contrast ? 'TRUE' : 'FALSE'
    def route = "--stage '${params.stage}' --process-fastq ${params.process_fastq ? 'TRUE' : 'FALSE'} --raw-samples '${params.raw_sample_names.join(',')}'"
    def orders = 'historical_column_order' in params.atac_modes ? "--header-order '${params.historical_header_samples.join(',')}' --data-order '${params.historical_data_samples.join(',')}'" : ''
    """
    validate_inputs.py --samplesheet "${samplesheet}" --source2target "${source2target}" \
      --reverse-contrast ${reverse} ${orders} ${route} --output validated_inputs.json
    """
}
