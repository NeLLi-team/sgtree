process WRITE_COLOR_FILE {
    publishDir params.outdir, mode: 'copy'

    input:
    val genomedir
    val refdir

    output:
    path 'color.txt', emit: color

    script:
    def ref_arg = refdir ? "--refdir '${refdir}'" : ''
    """
    write_color_file.py --genomedir '${genomedir}' ${ref_arg} --out color.txt
    """
}

process RENDER_TREE {
    publishDir params.outdir, mode: 'copy'

    input:
    path tree
    path color_file
    val  out_png

    output:
    path "${out_png}", emit: png

    script:
    """
    render_tree.py --tree ${tree} --color ${color_file} --out ${out_png}
    """
}
