process PCGR {
    tag "${meta.patient}:${meta.sample}"
    label 'process_high'

    conda "pcgr::pcgr=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:2.2.1':
        'docker.io/sigven/pcgr:2.2.1' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"
    path pon
    path vep_cache 

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = './PCGR'
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def cna      = params.cna_analysis ? "--input_cna $cna" : ''
    """
    export XDG_CACHE_HOME=/tmp  
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p $prefix

    pcgr \\
        --input_vcf $vcf \\
        --vep_dir $vep_cache \\
        --refdata_dir $database \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $prefix \\
        --tumor_dp_tag 'TDP' \\
        --tumor_af_tag 'TAF' \\
        --control_dp_tag 'NDP' \\
        --control_af_tag 'NAF' \\
        --call_conf_tag 'TAL' \\
        $cna \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr //g' ))
    END_VERSIONS
    """
}
