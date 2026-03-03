include { SGTREE_MAIN } from './sgtree_main'

workflow PREPARE_REFERENCE {
    take:
    refdir
    modeldir
    percent_models

    main:
    // Run the core pipeline on reference genomes (no marker selection, no ref-of-ref)
    SGTREE_MAIN(
        refdir,
        modeldir,
        percent_models,
        0,            // lflt
        'hmmalign',   // aln method
        false,        // has_ref
        file('NO_REF_MERGED'),
        file('NO_REF_PROTEOMES')
    )

    emit:
    table_elim_dups = SGTREE_MAIN.out.table_elim_dups
    proteomes       = SGTREE_MAIN.out.proteomes
}
