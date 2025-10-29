/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PCGR_GETREF            } from '../../modules/local/pcgr/getref'
include { ENSEMBLVEP_DOWNLOAD    } from '../../modules/nf-core/ensemblvep/download'

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
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTPRIORITIZATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input


    main:

    fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    if (params.database) { ch_pcgr_dir = Channel.fromPath("${params.database}/data/${params.genome.toLowerCase()}") } else { exit 1, "Please provide a path to the PCGR annotation database." }
    vep_cache               = Channel.fromPath(params.vep_cache)

    pcgr_bundle_version = params.pcgr_bundle_version ?: '20250314'

    //
    // GET REFERENCE DATA FOR PCGR
    //
    PCGR_GETREF([[id:'pcgr_reference'], pcgr_bundle_version, params.genome.toLowerCase()] )

    ch_pcgr_dir = PCGR_GETREF.out.pcgrref.map { _meta, pcgrref -> pcgrref }

    ENSEMBLVEP_DOWNLOAD([[id:'vep_cache'], params.genome, 'homo_sapiens', '113'])

    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map { _meta, cache -> cache}


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Preprocess VCF files
    //
    VCF_PREPROCESSING (
        ch_samplesheet,
        fasta
    )

    VCF_PREPROCESSING.out.filtered_ch
        .set { vcf_files }

    VCF_PREPROCESSING.out.ch_cna_files
        .set { cna_files }

    //
    // SUBWORKFLOW: Format input files
    //

    //vcf_files.view()
    //cna_files.view()

    FORMAT_FILES (
        vcf_files,
        cna_files
    )

    //
    // SUBWORKFLOW: pcgr
    //

    FORMAT_FILES.out.pcgr_ready_vcf.view()

    RUN_PCGR (
        FORMAT_FILES.out.pcgr_ready_vcf,
        ch_pcgr_dir.collect(),
        //FORMAT_FILES.out.pon_vcf,
        vep_cache.collect()
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'variantprioritization_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
