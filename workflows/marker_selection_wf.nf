include { TRIMAL_SIMPLE }                                        from '../modules/trim'
include { TRIMAL_SIMPLE as TRIMAL_SIMPLE_ROUND1 }                 from '../modules/trim'
include { TRIMAL_SIMPLE as TRIMAL_SIMPLE_FINAL }                  from '../modules/trim'
include { FASTTREE_PER_MARKER }                                   from '../modules/phylogeny'
include { FASTTREE_INTERMEDIATE }                                 from '../modules/phylogeny'
include { RF_SELECTION; REMOVE_SINGLES; WRITE_CLEANED_ALIGNMENT; FORMAT_RF_VALUES } from '../modules/marker_selection'
include { RF_SELECTION as RF_SELECTION_ROUND2 }                   from '../modules/marker_selection'
include { WRITE_CLEANED_ALIGNMENT as WRITE_CLEANED_ALIGNMENT_ROUND1 } from '../modules/marker_selection'
include { WRITE_CLEANED_ALIGNMENT as WRITE_CLEANED_ALIGNMENT_ROUND2 } from '../modules/marker_selection'
include { BUILD_SUPERMATRIX as BUILD_SUPERMATRIX_FINAL }          from '../modules/supermatrix'
include { BUILD_SUPERMATRIX as BUILD_SUPERMATRIX_ROUND1 }         from '../modules/supermatrix'
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
    global_rounds = (params.selection_global_rounds ?: 1) as int

    if (global_rounds > 1) {
        // Round 1 duplicate cleanup without singleton pruning.
        trees_with_aln_round1 = RF_SELECTION.out.cleaned_tree.join(aligned)
        WRITE_CLEANED_ALIGNMENT_ROUND1(trees_with_aln_round1)
        TRIMAL_SIMPLE_ROUND1(WRITE_CLEANED_ALIGNMENT_ROUND1.out.cleaned)

        cleaned_round1 = TRIMAL_SIMPLE_ROUND1.out.trimmed
            .map { marker, file -> file }
            .collect()

        BUILD_SUPERMATRIX_ROUND1(cleaned_round1)
        FASTTREE_INTERMEDIATE(BUILD_SUPERMATRIX_ROUND1.out.supermatrix)

        // Round 2 duplicate cleanup against the rebuilt guide tree.
        RF_SELECTION_ROUND2(
            FASTTREE_PER_MARKER.out.tree,
            FASTTREE_INTERMEDIATE.out.tree,
            table_elim_dups,
            ref_list
        )

        rf_values_raw = RF_SELECTION_ROUND2.out.rf_values.collectFile(name: 'rf_values_raw.txt')
        FORMAT_RF_VALUES(rf_values_raw)

        if (do_singles) {
            REMOVE_SINGLES(RF_SELECTION_ROUND2.out.cleaned_tree, FASTTREE_INTERMEDIATE.out.tree, table_elim_dups)
            trees_for_cleaning = REMOVE_SINGLES.out.tree
        } else {
            trees_for_cleaning = RF_SELECTION_ROUND2.out.cleaned_tree
        }

        trees_with_aln = trees_for_cleaning.join(aligned)
        WRITE_CLEANED_ALIGNMENT_ROUND2(trees_with_aln)
        TRIMAL_SIMPLE_FINAL(WRITE_CLEANED_ALIGNMENT_ROUND2.out.cleaned)
    } else {
        rf_values_raw = RF_SELECTION.out.rf_values.collectFile(name: 'rf_values_raw.txt')
        FORMAT_RF_VALUES(rf_values_raw)

        if (do_singles) {
            REMOVE_SINGLES(RF_SELECTION.out.cleaned_tree, species_tree, table_elim_dups)
            trees_for_cleaning = REMOVE_SINGLES.out.tree
        } else {
            trees_for_cleaning = RF_SELECTION.out.cleaned_tree
        }

        trees_with_aln = trees_for_cleaning.join(aligned)
        WRITE_CLEANED_ALIGNMENT_ROUND1(trees_with_aln)
        TRIMAL_SIMPLE_FINAL(WRITE_CLEANED_ALIGNMENT_ROUND1.out.cleaned)
    }

    cleaned_collected = TRIMAL_SIMPLE_FINAL.out.trimmed
        .map { marker, file -> file }
        .collect()

    BUILD_SUPERMATRIX_FINAL(cleaned_collected)
    FASTTREE_FINAL(BUILD_SUPERMATRIX_FINAL.out.supermatrix)

    emit:
    tree_final = FASTTREE_FINAL.out.tree
    rf_values  = FORMAT_RF_VALUES.out.rf_values
}
