process CPSR {
    tag "${meta.patient}.${meta.sample}"
    label 'process_medium'

    conda "pcgr::pcgr=2.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sigven/pcgr:2.1.2':
        'docker.io/sigven/pcgr:2.1.2' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"

    output:
    tuple val(meta), path("${meta.patient}.${meta.sample}"), emit: cpsr_reports
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = "./PCGR"
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    mkdir -p $prefix

    cpsr \\
        --input_vcf $vcf \\
        --pcgr_dir $database \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( cpsr --version | sed 's/cpsr //g' ))
    END_VERSIONS
    """
}
