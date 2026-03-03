process CONCAT_INPUTS {
    input:
    path genomedir
    path modeldir

    output:
    path 'models.hmm'    , emit: models
    path 'models_split'  , emit: models_split
    path 'proteomes.faa' , emit: proteomes
    env  MODEL_COUNT     , emit: model_count

    script:
    """
    if [ -d "${modeldir}" ]; then
        cat ${modeldir}/*.hmm > models.hmm
    else
        cp ${modeldir} models.hmm
    fi
    if [ -d "${genomedir}" ]; then
        cat ${genomedir}/*.faa > proteomes.faa
    else
        cp ${genomedir} proteomes.faa
    fi
    split_hmm_models.py --models models.hmm --outdir models_split
    MODEL_COUNT=\$(grep -c '^NAME' models.hmm)
    """
}

process HMMSEARCH {
    cpus { params.hmmsearch_cpus ?: 8 }

    input:
    path models
    path proteomes

    output:
    path 'hits.hmmout', emit: hmmout

    script:
    """
    run_hmmsearch_pyhmmer.py \\
        --models ${models} \\
        --proteomes ${proteomes} \\
        --out hits.hmmout \\
        --cpus ${task.cpus} \\
        --hmmsearch_cutoff ${params.hmmsearch_cutoff ?: 'cut_ga'} \\
        --hmmsearch_evalue ${params.hmmsearch_evalue ?: 1e-5}
    """
}
