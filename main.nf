#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SGTREE_MAIN }         from './workflows/sgtree_main'
include { MARKER_SELECTION_WF } from './workflows/marker_selection_wf'
include { PREPARE_REFERENCE }   from './workflows/prepare_reference'
include { WRITE_HEATMAP }       from './modules/itol'
include { WRITE_COLOR_FILE; RENDER_TREE } from './modules/render'

workflow {
    // Validate required params
    if (!params.genomedir) error "ERROR: --genomedir is required"
    if (!params.modeldir)  error "ERROR: --modeldir is required"

    genomedir = file(params.genomedir, checkIfExists: true)
    modeldir  = file(params.modeldir, checkIfExists: true)
    color_refdir = params.ref ? file(params.ref, checkIfExists: true).toString() : ''
    marker_selection_enabled = (
        (params.marker_selection instanceof Boolean)
            ? params.marker_selection
            : ['true', 'yes', '1'].contains(params.marker_selection?.toString()?.trim()?.toLowerCase())
    )
    singles_enabled = (
        (params.singles instanceof Boolean)
            ? params.singles
            : ['true', 'yes', '1'].contains(params.singles?.toString()?.trim()?.toLowerCase())
    )

    // Legacy color annotation artifact (query=red, ref=gray)
    WRITE_COLOR_FILE(genomedir.toString(), color_refdir)

    // Handle reference directory
    if (params.ref) {
        refdir = file(params.ref, checkIfExists: true)

        // Run reference sub-workflow
        PREPARE_REFERENCE(
            refdir,
            modeldir,
            params.percent_models,
            params.max_sdup,
            params.max_dupl
        )

        ref_merged_final = PREPARE_REFERENCE.out.table_elim_dups
        ref_proteomes    = PREPARE_REFERENCE.out.proteomes

        // Build ref list for marker selection (filenames like IMG123.faa)
        ref_list = Channel.fromPath("${params.ref}/*.faa")
            .map { it.name }
            .collect()
            .map { it.join(',') }
    } else {
        ref_merged_final = file('NO_REF_MERGED')
        ref_proteomes    = file('NO_REF_PROTEOMES')
        ref_list         = Channel.value('')
    }

    // Run core pipeline (steps 1-9)
    SGTREE_MAIN(
        genomedir,
        modeldir,
        params.percent_models,
        params.lflt,
        params.max_sdup,
        params.max_dupl,
        params.aln,
        params.ref ? true : false,
        ref_merged_final,
        ref_proteomes
    )

    // Marker selection (steps 10-15)
    if (marker_selection_enabled) {
        MARKER_SELECTION_WF(
            SGTREE_MAIN.out.aligned,
            SGTREE_MAIN.out.tree,
            SGTREE_MAIN.out.table_elim_dups,
            ref_list,
            singles_enabled
        )

        // iTOL heatmap on final tree
        WRITE_HEATMAP(
            MARKER_SELECTION_WF.out.tree_final,
            SGTREE_MAIN.out.marker_count_matrix,
            'marker_counts.txt'
        )

        // Legacy final rendered tree artifact
        if (params.render_png) {
            RENDER_TREE(
                MARKER_SELECTION_WF.out.tree_final,
                WRITE_COLOR_FILE.out.color,
                'tree_final.png'
            )
        }
    } else {
        // iTOL heatmap on basic tree
        WRITE_HEATMAP(
            SGTREE_MAIN.out.tree,
            SGTREE_MAIN.out.marker_count_matrix,
            'marker_count.txt'
        )
    }
}
