process BUILD_SUPERMATRIX {
    input:
    path trimmed_files  // collected from all markers

    output:
    path 'concatenated.faa', emit: supermatrix
    path 'table_df_concatenated_w_X', emit: table

    script:
    """
    build_supermatrix.py --trimmed_dir . --out concatenated.faa --table table_df_concatenated_w_X
    """
}
