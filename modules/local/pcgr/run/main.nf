process PCGR_RUN {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13d45cfcdc2a705e5fc5e8eba1a27c54aecba3f483a085aa13505560a2528932/data'
        : 'community.wave.seqera.io/library/pcgr:2.2.5.9000--2d18ab4fab1bbfe0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path pcgr_dir
    //path pon
    path vep_cache

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    tuple val("${task.process}"), val('pcgr'),  eval("pcgr --version | sed 's/pcgr //g'"), topic: versions, emit: versions_pcgr

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome  = task.ext.genome ?: ''
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def cna_cmd = params.cna_analysis ? "--input_cna ${cna}" : ''
    """
    export XDG_CACHE_HOME=/tmp
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p ${prefix}

    pcgr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache} \\
        --refdata_dir ${pcgr_dir} \\
        --output_dir ${prefix} \\
        --genome_assembly ${genome} \\
        --sample_id ${prefix} \\
        --tumor_dp_tag 'TDP' \\
        --tumor_af_tag 'TAF' \\
        --call_conf_tag 'TAL' \\
        ${cna_cmd} \\
        ${args}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    """
}
