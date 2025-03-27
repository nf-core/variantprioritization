process PCGR_VCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/barryd237/pysam-xcmds:latest' :
        'docker.io/barryd237/pysam-xcmds:latest' }"

    input:
    tuple val(meta), path(keys), path(vcf), path(tbi)
    path(pcgr_header)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pcgr_vcf.py \\
        --sample ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$( python --version | cut -d' ' -f2 ))
    END_VERSIONS
    """
}
