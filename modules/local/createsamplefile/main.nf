process CREATESAMPLEFILE {
    tag "${meta.id}"
    label 'process_single'

    input:
    val(meta)

    output:
    tuple val(meta), path("${filename}"), emit: samplefile
    tuple val("${task.process}"), val('createsamplefile'),  eval("echo $VERSION"), topic: versions, emit: versions_createsamplefile

    script:
    def prefix = task.ext.prefix   ?: "${meta.id}"
    filename = task.ext.prefix ? "${prefix}.txt" : "${meta.id}.${meta.tool}.txt"
    VERSION = '1.0.0'
    """
    echo "${meta.id}" > ${filename}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    filename = task.ext.prefix ? "${prefix}.txt" : "${meta.id}.${meta.tool}.txt"
    """
    touch ${filename}
    """
}