process CONCAT_INPUTS {
    input:
    path genomedir
    path modeldir

    output:
    path 'models.hmm'    , emit: models
    path 'proteomes.faa' , emit: proteomes
    env  MODEL_COUNT     , emit: model_count

    script:
    """
    cat ${modeldir}/*.hmm > models.hmm
    if [ -d "${genomedir}" ]; then
        cat ${genomedir}/*.faa > proteomes.faa
    else
        cp ${genomedir} proteomes.faa
    fi
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
    hmmsearch --cut_ga --cpu ${task.cpus} --domtblout hits.hmmout --noali ${models} ${proteomes} > /dev/null
    """
}
