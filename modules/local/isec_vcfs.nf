process ISEC_SOMATIC_VCFS {
    tag "${meta.patient}:${meta.sample}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/barryd237/pysam-xcmds:latest' :
        'docker.io/barryd237/pysam-xcmds:latest' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}_keys.txt"), emit: variant_tool_map
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}" // meta.sample, toggle using modules.config
    """
    isec_vcfs.py \
        --sample ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$( python --version | cut -d' ' -f2 ))
    END_VERSIONS
    """
}
