process CPSR_VALIDATE_INPUT {
    tag "$meta.id"
    label 'process_medium'

    conda "pcgr::pcgr=2.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:2.1.2':
        'docker.io/sigven/pcgr:2.1.2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"

    output:
    tuple val(meta), path("${prefix}/*cpsr_ready.vcf"), emit: validated_vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def database = "./PCGR"
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    mkdir $prefix

    cpsr_validate_input.py \\
        $database \\
        $vcf \\
        "None" \\
        "None" \\
        0 \\
        ${params.genome.toLowerCase()} \\
        $prefix \\
        $params.panel_id \\
        0 \\
        --output_dir $prefix \\
        --debug

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr //g' ))
    END_VERSIONS
    """
}
