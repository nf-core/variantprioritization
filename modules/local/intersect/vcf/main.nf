process INTERSECT_VCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f7a55b8b6c7e72be975e8dfd6cba9117f9c6ec41e67e5c3a095e6021fad7d71/data'
        : 'community.wave.seqera.io/library/bcftools_pandas_python:00d964080a86c66f'}"

    input:
    tuple val(meta), path(isec_results)

    output:
    tuple val(meta), path("${prefix}_keys.txt"), emit: variant_tool_map
    tuple val("${task.process}"), val('python'),  eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // meta.sample, toggle using modules.config
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    isec_vcfs.py \
        --sample ${prefix}
    """
}
