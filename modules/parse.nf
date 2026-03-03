process PARSE_HMMSEARCH {
    publishDir params.outdir, mode: 'copy', pattern: 'marker_count_matrix.csv'
    publishDir params.outdir, mode: 'copy', pattern: 'log_genomes_removed.txt'

    input:
    path hmmout
    path proteomes, stageAs: 'query_proteomes.faa'
    val  model_count
    val  percent_models
    val  lflt
    val  max_sdup
    val  max_dupl
    val  has_ref
    path ref_merged_final, stageAs: 'ref_merged_final.csv'
    path ref_proteomes, stageAs: 'ref_proteomes.faa'

    output:
    path 'extracted/*'       , emit: marker_id_lists
    path 'table_elim_dups'   , emit: table_elim_dups
    path 'combined_proteomes.faa', emit: combined_proteomes
    path 'combined_proteomes.idx', emit: combined_proteomes_idx
    path 'marker_count_matrix.csv', emit: marker_count_matrix
    path 'log_genomes_removed.txt', emit: log_removed

    script:
    def ref_merged_arg = has_ref ? "--ref_merged_final ${ref_merged_final}" : ''
    def ref_prot_arg   = has_ref ? "--ref_proteomes ${ref_proteomes}" : ''
    """
    parse_hmmsearch.py \\
        --hmmout ${hmmout} \\
        --proteomes ${proteomes} \\
        --model_count ${model_count} \\
        --percent_models ${percent_models} \\
        --lflt ${lflt} \\
        --max_sdup ${max_sdup} \\
        --max_dupl ${max_dupl} \\
        ${ref_merged_arg} \\
        ${ref_prot_arg}
    """
}
