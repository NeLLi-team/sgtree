process TRIMAL {
    tag "${marker}"

    input:
    tuple val(marker), path(alignment)

    output:
    tuple val(marker), path("trimmed_${marker}.faa"), emit: trimmed

    script:
    """
    trimal -in ${alignment} -out trimmed_${marker}.faa -gt 0.1
    """
}

process TRIMAL_SIMPLE {
    tag "${marker}"

    input:
    tuple val(marker), path(alignment)

    output:
    tuple val(marker), path("trimmed_${marker}.faa"), emit: trimmed

    script:
    """
    trimal -in ${alignment} -out trimmed_${marker}.faa -gt 0.1
    """
}
