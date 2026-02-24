process CREATESAMPLEFILE {
    tag "${meta.id}"
    label 'process_single'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: samplefile
    tuple val("${task.process}"), val('createsamplefile'), eval("echo ${VERSION}"), topic: versions, emit: versions_createsamplefile

    script:
    VERSION = '1.0.0'
    """
    echo "${meta.sample}" > ${meta.id}.txt
    """

    stub:
    """
    touch ${meta.id}.txt
    """
}
