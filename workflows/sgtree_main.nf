include { CONCAT_INPUTS; HMMSEARCH } from '../modules/search'
include { ANI_CLUSTER }              from '../modules/ani'
include { PARSE_HMMSEARCH }          from '../modules/parse'
include { EXTRACT_SEQUENCES }        from '../modules/extract'
include { ALIGN_MAFFT; ALIGN_MAFFT_LINSI; ALIGN_HMMALIGN } from '../modules/align'
include { ELIMINATE_DUPLICATES }     from '../modules/duplicates'
include { TRIMAL }                   from '../modules/trim'
include { BUILD_SUPERMATRIX }        from '../modules/supermatrix'
include { FASTTREE }                 from '../modules/phylogeny'

workflow SGTREE_MAIN {
    take:
    genomedir
    modeldir
    percent_models
    lflt
    max_sdup
    max_dupl
    aln_method
    ani_enabled
    has_ref
    ref_merged_final   // path or NO_REF placeholder
    ref_proteomes      // path or NO_REF placeholder
    ref_genome_manifest // path or NO_REF placeholder

    main:
    // Step 1: Concatenate inputs
    CONCAT_INPUTS(genomedir, modeldir)

    if (ani_enabled) {
        ANI_CLUSTER(
            CONCAT_INPUTS.out.genome_manifest,
            has_ref,
            ref_genome_manifest
        )
        keep_genomes = ANI_CLUSTER.out.kept
        ani_clusters_out = ANI_CLUSTER.out.clusters
        ani_representatives_out = ANI_CLUSTER.out.representatives
        ani_pairs_out = ANI_CLUSTER.out.pairs
    } else {
        keep_genomes = file('NO_KEEP_GENOMES')
        ani_clusters_out = file('NO_ANI_CLUSTERS')
        ani_representatives_out = file('NO_ANI_REPRESENTATIVES')
        ani_pairs_out = file('NO_ANI_PAIRS')
    }

    // Step 2: Run hmmsearch
    HMMSEARCH(CONCAT_INPUTS.out.models, CONCAT_INPUTS.out.proteomes)

    // Step 3-4: Parse results, extract per-marker ID lists, build combined proteomes
    PARSE_HMMSEARCH(
        HMMSEARCH.out.hmmout,
        CONCAT_INPUTS.out.proteomes,
        CONCAT_INPUTS.out.model_count,
        percent_models,
        lflt,
        max_sdup,
        max_dupl,
        ani_enabled,
        keep_genomes,
        has_ref,
        ref_merged_final,
        ref_proteomes
    )

    // Scatter per marker: create [marker, id_list] tuples
    marker_ch = PARSE_HMMSEARCH.out.marker_id_lists
        .flatten()
        .map { file -> [ file.baseName, file ] }

    // Step 4: Extract sequences per marker
    EXTRACT_SEQUENCES(
        marker_ch,
        PARSE_HMMSEARCH.out.combined_proteomes,
        PARSE_HMMSEARCH.out.combined_proteomes_idx
    )

    // Step 5: Alignment (method-dependent)
    if (aln_method == 'mafft') {
        ALIGN_MAFFT(EXTRACT_SEQUENCES.out.seqs)
        aligned = ALIGN_MAFFT.out.aligned
    } else if (aln_method == 'mafft-linsi') {
        ALIGN_MAFFT_LINSI(EXTRACT_SEQUENCES.out.seqs)
        aligned = ALIGN_MAFFT_LINSI.out.aligned
    } else {
        ALIGN_HMMALIGN(EXTRACT_SEQUENCES.out.seqs, CONCAT_INPUTS.out.models_split)
        aligned = ALIGN_HMMALIGN.out.aligned
    }

    // Step 6: Eliminate duplicates
    ELIMINATE_DUPLICATES(aligned, PARSE_HMMSEARCH.out.table_elim_dups)

    // Step 7: Trim alignments
    TRIMAL(ELIMINATE_DUPLICATES.out.deduped)

    // Step 8: Build supermatrix (collect all trimmed files)
    trimmed_collected = TRIMAL.out.trimmed.map { marker, file -> file }.collect()
    BUILD_SUPERMATRIX(trimmed_collected)

    // Step 9: Build species tree
    FASTTREE(BUILD_SUPERMATRIX.out.supermatrix)

    emit:
    tree               = FASTTREE.out.tree
    aligned            = aligned
    table_elim_dups    = PARSE_HMMSEARCH.out.table_elim_dups
    marker_count_matrix = PARSE_HMMSEARCH.out.marker_count_matrix
    proteomes          = CONCAT_INPUTS.out.proteomes
    genome_manifest    = CONCAT_INPUTS.out.genome_manifest
    ani_clusters       = ani_clusters_out
    ani_representatives = ani_representatives_out
    ani_pairs          = ani_pairs_out
}
