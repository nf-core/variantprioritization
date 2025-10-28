process PREPARE_PCGRREF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'community.wave.seqera.io/library/coreutils_grep_gzip_tar_wget:5509a0d21b41d5be' }"

    input:
    tuple val(meta), val(bundleversion), val(genome)

    output:
    tuple val(meta), path("${bundleversion}"), emit: pcgrref
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bundle = "pcgr_ref_data.${bundleversion}.${genome}.tgz"
    """
    # BUNDLE_VERSION="20250314"
    # GENOME="grch38" # or "grch37"
    wget https://insilico.hpc.uio.no/pcgr/${bundle}
    gzip -dc ${bundle} | tar xvf -

    mkdir ${bundleversion}
    mv data/ ${bundleversion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare: \$(prepare --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir ${bundleversion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare: \$(prepare --version)
    END_VERSIONS
    """
}
