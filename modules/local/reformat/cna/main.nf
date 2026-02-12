process REFORMAT_CNA {
    tag "${meta.patient}:${meta.sample}:${cna}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/935400e4df07528697155f191c225dbf18ac4dc5d7779b1b14f5b974e9237227/data'
        : 'community.wave.seqera.io/library/pysam_vcf2tsvpy_numpy_pandas_python:eb0ee661861e1b56'}"

    input:
    tuple val(meta), path(cna)

    output:
    tuple val(meta), path("${prefix}.*.tsv"), emit: cna
    tuple val("${task.process}"), val('python'),  eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.sample}"
    """
    reformat_cna.py \\
        --input ${cna} \\
        --sample ${prefix}
    """
}
