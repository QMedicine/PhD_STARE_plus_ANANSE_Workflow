process HIC_CONTACTS_FROM_H5 {
    tag state
    label 'hicexplorer'
    cpus 4
    memory '16 GB'
    time '6h'
    input:
    tuple val(state), path(matrix)
    output:
    tuple val(state), path('contacts'), emit: contacts
    script:
    """
    hicConvertFormat --matrices "${matrix}" --inputFormat h5 --outputFormat ginteractions --outFileName contacts_all
    split_contacts.py --input contacts_all.tsv --outdir contacts
    """
}

process STARE_GABC {
    tag "${mode}:${state}"
    label 'stare'
    cpus 32
    memory '128 GB'
    time '2h'
    publishDir { "${params.outdir}/${mode}/stare/${state}" }, mode: 'copy'
    input:
    tuple val(mode), val(state), val(columns), path(activity), path(contacts)
    path samplesheet
    path fasta
    path gtf
    path psem
    path modified_stare
    output:
    tuple val(mode), val(state), path("${state}/Gene_TF_matrices/*TF_Gene_Affinities*.txt.gz"), emit: affinities
    path "${state}", emit: results
    script:
    // Copy the exact shell script alongside symlinks to the container's compiled binaries.
    // This preserves its dirname($0) executable lookup without changing the archived script.
    """
    mkdir stare_bin
    cp "${modified_stare}" stare_bin/STARE.sh
    for executable in ReplaceInvalidChars TRAPmulti STARE_ABCpp TF_Gene_Scorer_Reshape PSCM_to_PSEM Promoter_Windows; do
      ln -s "/usr/local/bin/\$executable" "stare_bin/\$executable"
    done
    OMP_NUM_THREADS=${task.cpus} bash stare_bin/STARE.sh -b "${activity}" -n "${columns}" -a "${gtf}" -g "${fasta}" \
      -f "${contacts}" -k ${params.chosen_res} -w ${params.stare_window} -m ${params.enhancer_window} \
      -t ${params.gabc_threshold} -p "${psem}" -c ${task.cpus} -i all_tss -q True -d True -o '${state}'
    """
}
