include { TRIMAL_SIMPLE }                                        from '../modules/trim'
include { TRIMAL_SIMPLE as TRIMAL_SIMPLE_FINAL }                  from '../modules/trim'
include { FASTTREE_PER_MARKER }                                   from '../modules/phylogeny'
include { RF_SELECTION; REMOVE_SINGLES; WRITE_CLEANED_ALIGNMENT; FORMAT_RF_VALUES } from '../modules/marker_selection'
include { BUILD_SUPERMATRIX as BUILD_SUPERMATRIX_FINAL }          from '../modules/supermatrix'
include { FASTTREE_FINAL }                                        from '../modules/phylogeny'

workflow MARKER_SELECTION_WF {
    take:
    aligned            // channel: [marker, aligned.faa]
    species_tree       // path: tree.nwk
    table_elim_dups    // path: table_elim_dups
    ref_list           // val: comma-separated ref filenames or empty string
    do_singles         // val: boolean

    main:
    // Step 10: Trim aligned files (with duplicates) for protein trees
    TRIMAL_SIMPLE(aligned)

    // Step 11: Build per-marker protein trees
    FASTTREE_PER_MARKER(TRIMAL_SIMPLE.out.trimmed)

    // Step 12: RF-distance marker selection
    RF_SELECTION(
        FASTTREE_PER_MARKER.out.tree,
        species_tree,
        table_elim_dups,
        ref_list
    )

    // Collect per-marker RF values and add legacy header line
    rf_values_raw = RF_SELECTION.out.rf_values.collectFile(name: 'rf_values_raw.txt')
    FORMAT_RF_VALUES(rf_values_raw)

    // Step 13: Optionally remove singletons
    if (do_singles) {
        REMOVE_SINGLES(RF_SELECTION.out.cleaned_tree, species_tree)
        trees_for_cleaning = REMOVE_SINGLES.out.tree
    } else {
        trees_for_cleaning = RF_SELECTION.out.cleaned_tree
    }

    // Step 14: Write cleaned alignments (join tree back with original alignment by marker)
    trees_with_aln = trees_for_cleaning.join(aligned)
    WRITE_CLEANED_ALIGNMENT(trees_with_aln)

    // Step 15: Trim cleaned alignments, build final supermatrix and tree
    TRIMAL_SIMPLE_FINAL(WRITE_CLEANED_ALIGNMENT.out.cleaned)

    cleaned_collected = TRIMAL_SIMPLE_FINAL.out.trimmed
        .map { marker, file -> file }
        .collect()

    BUILD_SUPERMATRIX_FINAL(cleaned_collected)
    FASTTREE_FINAL(BUILD_SUPERMATRIX_FINAL.out.supermatrix)

    emit:
    tree_final = FASTTREE_FINAL.out.tree
    rf_values  = FORMAT_RF_VALUES.out.rf_values
}
