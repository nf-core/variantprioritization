/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFERENCE_DATA         } from '../../subworkflows/local/reference_data'

include { VCF_PREPROCESSING      } from '../../subworkflows/local/vcf_preprocessing'
include { FORMAT_FILES           } from '../../subworkflows/local/format_files'
include { PCGR as RUN_PCGR       } from '../../modules/local/pcgr'

include { getGenomeAttribute     } from '../../subworkflows/local/utils_nfcore_variantprioritization_pipeline'
include { MULTIQC                } from '../../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_variantprioritization_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTPRIORITIZATION {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    params

    main:
    params.fasta = getGenomeAttribute('fasta')
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def fasta = params.fasta ? channel.fromPath(params.fasta, checkIfExists: true).map { it -> [[id: it.baseName], it] }.collect() : channel.empty()

    REFERENCE_DATA(
        params.genome,
        params.pcgr_download,
        params.pcgr_bundleversion,
        params.pcgr_database_dir,
        params.vep_cache,
        params.vep_cache_version,
        params.vep_species,
    )

    def ch_pcgr_dir = REFERENCE_DATA.out.pcgr_dir
    def ch_vep_cache = REFERENCE_DATA.out.vep_cache

    ch_versions = ch_versions.mix(REFERENCE_DATA.out.versions)

    //
    // SUBWORKFLOW: Preprocess VCF files
    //
    VCF_PREPROCESSING(
        ch_samplesheet,
        fasta,
    )

    def ch_vcf_files = VCF_PREPROCESSING.out.filtered_ch
    def ch_cna_files = VCF_PREPROCESSING.out.ch_cna_files

    ch_versions = ch_versions.mix(VCF_PREPROCESSING.out.versions)

    //
    // SUBWORKFLOW: Format input files
    //
    FORMAT_FILES(
        ch_vcf_files,
        ch_cna_files,
    )

    ch_versions = ch_versions.mix(FORMAT_FILES.out.versions)

    //
    // SUBWORKFLOW: pcgr
    //
    RUN_PCGR(
        FORMAT_FILES.out.pcgr_ready_vcf,
        ch_pcgr_dir.collect(),
        ch_vep_cache.collect(),
    )

    ch_versions = ch_versions.mix(RUN_PCGR.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'variantprioritization_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
