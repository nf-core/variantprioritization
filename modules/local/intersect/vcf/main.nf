process INTERSECT_VCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e865b57ba6a9b7164c8018cf631df0ae2746cf5ca3db5666502fc0d61d9bcf91/data'
        : 'community.wave.seqera.io/library/bcftools_pysam_pandas_python:6b813c53a7ef4ede'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}_keys.txt"), emit: variant_tool_map
    tuple val("${task.process}"), val('python'),  eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    // meta.sample, toggle using modules.config
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    isec_vcfs.py \\
        --sample ${prefix}
    """
}
