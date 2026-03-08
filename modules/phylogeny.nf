process FASTTREE {
    publishDir params.outdir, mode: 'copy'
    cpus { params.fasttree_cpus ?: 1 }

    input:
    path supermatrix

    output:
    path 'tree.nwk', emit: tree

    script:
    """
    method="${params.tree_method ?: 'fasttree'}"
    iqfast="${params.iqtree_fast ?: true}"
    iqmodel="${params.iqtree_model ?: 'LG+F+I+G4'}"
    iqfast_norm=\$(echo "\${iqfast}" | tr '[:upper:]' '[:lower:]')
    if [ "\${method}" = "iqtree" ]; then
        iq_extra=""
        if [ "\${iqfast_norm}" = "true" ] || [ "\${iqfast_norm}" = "yes" ] || [ "\${iqfast_norm}" = "1" ]; then
            iq_extra="-fast"
        fi
        iqtree --quiet --prefix tree_iq -m "\${iqmodel}" -T ${task.cpus} \${iq_extra} -s ${supermatrix}
        cp tree_iq.treefile tree.nwk
    else
        FastTree -quiet -out tree.nwk ${supermatrix}
    fi
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
    method="${params.tree_method ?: 'fasttree'}"
    iqfast="${params.iqtree_fast ?: true}"
    iqmodel="${params.iqtree_model ?: 'LG+F+I+G4'}"
    iqfast_norm=\$(echo "\${iqfast}" | tr '[:upper:]' '[:lower:]')
    if [ "\${method}" = "iqtree" ]; then
        iq_extra=""
        if [ "\${iqfast_norm}" = "true" ] || [ "\${iqfast_norm}" = "yes" ] || [ "\${iqfast_norm}" = "1" ]; then
            iq_extra="-fast"
        fi
        iqtree --quiet --prefix tree_final_iq -m "\${iqmodel}" -T ${task.cpus} \${iq_extra} -s ${supermatrix}
        cp tree_final_iq.treefile tree_final.nwk
    else
        FastTree -quiet -out tree_final.nwk ${supermatrix}
    fi
    """
}

process FASTTREE_INTERMEDIATE {
    cpus { params.fasttree_cpus ?: 1 }

    input:
    path supermatrix

    output:
    path 'tree_round.nwk', emit: tree

    script:
    """
    method="${params.tree_method ?: 'fasttree'}"
    iqfast="${params.iqtree_fast ?: true}"
    iqmodel="${params.iqtree_model ?: 'LG+F+I+G4'}"
    iqfast_norm=\$(echo "\${iqfast}" | tr '[:upper:]' '[:lower:]')
    if [ "\${method}" = "iqtree" ]; then
        iq_extra=""
        if [ "\${iqfast_norm}" = "true" ] || [ "\${iqfast_norm}" = "yes" ] || [ "\${iqfast_norm}" = "1" ]; then
            iq_extra="-fast"
        fi
        iqtree --quiet --prefix tree_round_iq -m "\${iqmodel}" -T ${task.cpus} \${iq_extra} -s ${supermatrix}
        cp tree_round_iq.treefile tree_round.nwk
    else
        FastTree -quiet -out tree_round.nwk ${supermatrix}
    fi
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
    method="${params.tree_method ?: 'fasttree'}"
    iqfast="${params.iqtree_fast ?: true}"
    iqmodel="${params.iqtree_model ?: 'LG+F+I+G4'}"
    iqfast_norm=\$(echo "\${iqfast}" | tr '[:upper:]' '[:lower:]')
    if [ "\${method}" = "iqtree" ]; then
        iq_extra=""
        if [ "\${iqfast_norm}" = "true" ] || [ "\${iqfast_norm}" = "yes" ] || [ "\${iqfast_norm}" = "1" ]; then
            iq_extra="-fast"
        fi
        iqtree --quiet --prefix marker_iq -m "\${iqmodel}" -T ${task.cpus} \${iq_extra} -s ${alignment}
        cp marker_iq.treefile ${marker}_tree.out
    else
        FastTree -quiet -out ${marker}_tree.out ${alignment}
    fi
    """
}
