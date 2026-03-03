process EXTRACT_SEQUENCES {
    tag "${marker}"

    input:
    tuple val(marker), path(id_list)
    path proteomes

    output:
    tuple val(marker), path("${marker}.faa"), emit: seqs

    script:
    """
    extract_sequences.py --id_list ${id_list} --proteomes ${proteomes} --out ${marker}.faa
    """
}
