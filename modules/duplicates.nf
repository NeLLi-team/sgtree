process ELIMINATE_DUPLICATES {
    tag "${marker}"

    input:
    tuple val(marker), path(alignment)
    path table_elim_dups

    output:
    tuple val(marker), path("deduped_${marker}.faa"), emit: deduped

    script:
    """
    eliminate_duplicates.py --alignment ${alignment} --table_elim_dups ${table_elim_dups} --out deduped_${marker}.faa
    """
}
