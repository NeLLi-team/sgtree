process ALIGN_MAFFT {
    tag "${marker}"
    cpus { params.align_cpus ?: 4 }

    input:
    tuple val(marker), path(seqs)

    output:
    tuple val(marker), path("${marker}.faa"), emit: aligned

    script:
    """
    mafft --auto --thread ${task.cpus} --quiet ${seqs} > ${marker}.faa
    """
}

process ALIGN_MAFFT_LINSI {
    tag "${marker}"
    cpus { params.align_cpus ?: 4 }

    input:
    tuple val(marker), path(seqs)

    output:
    tuple val(marker), path("${marker}.faa"), emit: aligned

    script:
    """
    mafft-linsi --thread ${task.cpus} --quiet ${seqs} > ${marker}.faa
    """
}

process ALIGN_HMMALIGN {
    tag "${marker}"

    input:
    tuple val(marker), path(seqs)
    path modelset

    output:
    tuple val(marker), path("${marker}.faa"), emit: aligned

    script:
    """
    run_hmmalign_pyhmmer.py --models ${modelset} --marker ${marker} --seqs ${seqs} --out ${marker}.faa --cpus ${task.cpus}
    """
}
