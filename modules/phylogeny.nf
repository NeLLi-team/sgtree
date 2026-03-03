process FASTTREE {
    publishDir params.outdir, mode: 'copy'
    cpus { params.fasttree_cpus ?: 1 }

    input:
    path supermatrix

    output:
    path 'tree.nwk', emit: tree

    script:
    """
    FastTree -quiet -out tree.nwk ${supermatrix}
    """
}

process FASTTREE_FINAL {
    publishDir params.outdir, mode: 'copy'
    cpus { params.fasttree_cpus ?: 1 }

    input:
    path supermatrix

    output:
    path 'tree_final.nwk', emit: tree

    script:
    """
    FastTree -quiet -out tree_final.nwk ${supermatrix}
    """
}

process FASTTREE_PER_MARKER {
    tag "${marker}"
    cpus { params.fasttree_cpus ?: 1 }

    input:
    tuple val(marker), path(alignment)

    output:
    tuple val(marker), path("${marker}_tree.out"), emit: tree

    script:
    """
    FastTree -quiet -out ${marker}_tree.out ${alignment}
    """
}
