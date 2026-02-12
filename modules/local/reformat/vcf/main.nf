process REFORMAT_VCF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/78186c6fc95e23c235ffba10839ffd751c7fd2a131bbd7d34d4ee4ec80edd784/data'
        : 'community.wave.seqera.io/library/bcftools_pysam_python:c6d15d978dc52fc5'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('python'),   eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    reformat_vcf.py \\
        --input ${vcf} \\
        --output ${prefix}.vcf
    """
}
