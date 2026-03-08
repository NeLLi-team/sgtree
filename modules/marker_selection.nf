process RF_SELECTION {
    tag "${marker}"

    input:
    tuple val(marker), path(marker_tree)
    path species_tree
    path table_elim_dups
    val  ref_list

    output:
    tuple val(marker), path("cleaned_${marker}.nw"), emit: cleaned_tree
    path 'rf_values.txt', emit: rf_values

    script:
    def refs_arg = ref_list ? "--ref_list '${ref_list}'" : ''
    """
    rf_marker_selection.py \\
        --marker_tree ${marker_tree} \\
        --species_tree ${species_tree} \\
        --table_elim_dups ${table_elim_dups} \\
        --marker ${marker} \\
        --out cleaned_${marker}.nw \\
        --rf_out rf_values.txt \\
        --selection_mode ${params.selection_mode ?: 'coordinate'} \\
        --selection_max_rounds ${params.selection_max_rounds ?: 5} \\
        --lock_references ${params.lock_references ?: false} \\
        ${refs_arg}
    """
}

process REMOVE_SINGLES {
    tag "${marker}"

    input:
    tuple val(marker), path(cleaned_tree)
    path species_tree
    path table_elim_dups

    output:
    tuple val(marker), path("no_singles_${marker}.nw"), emit: tree

    script:
    """
    remove_singles.py --tree ${cleaned_tree} --species_tree ${species_tree} --table_elim_dups ${table_elim_dups} --mode ${params.singles_mode ?: 'neighbor'} --num_nei ${params.num_nei ?: 0} --singles_min_rfdist ${params.singles_min_rfdist ?: 0.25} --out no_singles_${marker}.nw
    """
}

process WRITE_CLEANED_ALIGNMENT {
    tag "${marker}"

    input:
    tuple val(marker), path(newick), path(alignment)

    output:
    tuple val(marker), path("cleaned_${marker}.faa"), emit: cleaned

    script:
    """
    write_cleaned_alignment.py --newick ${newick} --alignment ${alignment} --out cleaned_${marker}.faa
    """
}

process FORMAT_RF_VALUES {
    publishDir params.outdir, mode: 'copy'

    input:
    path rf_values_raw

    output:
    path 'marker_selection_rf_values.txt', emit: rf_values

    script:
    """
    {
      echo "ProteinID MarkerGene RFdistance Status"
      cat ${rf_values_raw}
    } > marker_selection_rf_values.txt
    """
}
