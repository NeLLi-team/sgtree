process WRITE_HEATMAP {
    publishDir params.outdir, mode: 'copy'

    input:
    path tree
    path marker_count_matrix
    val  outsuffix

    output:
    path "${outsuffix}", emit: heatmap

    script:
    """
    write_itol_heatmap.py --tree ${tree} --matrix ${marker_count_matrix} --out ${outsuffix}
    """
}
