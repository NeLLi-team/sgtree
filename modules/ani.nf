process ANI_CLUSTER {
    publishDir params.outdir, mode: 'copy', pattern: 'ani/*'
    cpus { params.ani_cpus ?: 4 }

    input:
    path query_manifest
    val  has_ref
    path ref_manifest, stageAs: 'ref_manifest.tsv'

    output:
    path 'ani/ani_kept_genomes.txt', emit: kept
    path 'ani/ani_clusters.tsv', emit: clusters
    path 'ani/ani_representatives.tsv', emit: representatives
    path 'ani/ani_pairwise.tsv', emit: pairs

    script:
    def ref_arg = has_ref ? "--ref-manifest ${ref_manifest}" : ''
    """
    mkdir -p ani
    python ${projectDir}/bin/run_ani_clustering.py \\
        --query-manifest ${query_manifest} \\
        ${ref_arg} \\
        --outdir ani \\
        --ani-threshold ${params.ani_threshold ?: 95.0} \\
        --ani-backend ${params.ani_backend ?: 'auto'} \\
        --ani-mcl-inflation ${params.ani_mcl_inflation ?: 2.0} \\
        --num-cpus ${task.cpus}
    """
}

process BUILD_SNP_TREES {
    publishDir params.outdir, mode: 'copy', pattern: 'snp_trees/*'
    cpus { params.snp_cpus ?: params.ani_cpus ?: 4 }

    input:
    path clusters
    path query_manifest
    val  has_ref
    path ref_manifest, stageAs: 'ref_manifest.tsv'

    output:
    path 'snp_trees', emit: snp_trees

    script:
    def ref_arg = has_ref ? "--ref-manifest ${ref_manifest}" : ''
    """
    mkdir -p snp_trees
    python ${projectDir}/bin/build_snp_trees.py \\
        --clusters ${clusters} \\
        --query-manifest ${query_manifest} \\
        ${ref_arg} \\
        --outdir snp_trees \\
        --min-cluster-size ${params.snp_tree_min_cluster_size ?: 3} \\
        --num-cpus ${task.cpus}
    """
}
