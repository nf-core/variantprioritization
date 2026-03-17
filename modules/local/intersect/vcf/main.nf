process INTERSECT_VCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e865b57ba6a9b7164c8018cf631df0ae2746cf5ca3db5666502fc0d61d9bcf91/data'
        : 'community.wave.seqera.io/library/bcftools_pysam_pandas_python:6b813c53a7ef4ede'}"

    input:
    tuple val(meta), path(vcfs, stageAs: 'inputs/*'), path(tbis, stageAs: 'inputs/*')

    output:
    tuple val(meta), path("${prefix}_keys.txt"), emit: variant_tool_map
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam
    tuple val("${task.process}"), val('python'), eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def vcfList = vcfs instanceof List ? vcfs : [vcfs]
    def toolList = meta.tools instanceof List ? meta.tools : [meta.tools]
    def vcfFiles = vcfList.collect { vcfFile -> new File(vcfFile.toString()).name }.join(',')
    def toolNames = toolList.collect { toolName -> toolName.toString() }.join(',')
    """
    isec_vcfs.py \\
        --input_dir inputs \\
        --sample_files ${vcfFiles} \\
        --tool_names ${toolNames} \\
        --output ${prefix}_keys.txt
    """
}
