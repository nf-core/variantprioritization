process PCGR_PREPAREVCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/935400e4df07528697155f191c225dbf18ac4dc5d7779b1b14f5b974e9237227/data'
        : 'community.wave.seqera.io/library/pysam_vcf2tsvpy_numpy_pandas_python:eb0ee661861e1b56'}"

    input:
    tuple val(meta), path(keys), path(vcf), path(tbi)
    path pcgr_header

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('python'),  eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pcgr_vcf.py \\
        --sample ${prefix}
    """
}
