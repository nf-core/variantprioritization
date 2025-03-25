/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//This optionally tabix and bgzip the input VCF file [input_check]
include { TABIX_TABIX as TABIX_INPUT_VCF      } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIP_INPUT_VCF } from '../../modules/nf-core/tabix/bgziptabix/main'

//This reformat the input VCF file [format_files]
include { BCFTOOLS_NORM as NORMALISE_VARIANTS } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER as FILTER_VARIANTS  } from '../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX as TABIX_FILTERED       } from '../../modules/nf-core/tabix/tabix/main'
//include { REFORMAT_VCF                      } from '../../modules/local/pcgr_reformat'
//include { REFORMAT_CNA                      } from '../../modules/local/pcgr_reformat'
//include { REFORMAT_PON                      } from '../../modules/local/pcgr_reformat' 

//This merge the files [merge_vcfs]
include { TABIX_BGZIPTABIX as BGZIPTABIX_CPSR               } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_CONCAT                                   } from '../../modules/nf-core/bcftools/concat/main'
include { TABIX_TABIX as TABIX_CONCAT                       } from '../../modules/nf-core/tabix/tabix/main'
//include { ISEC_SOMATIC_VCFS as INTERSECT_SOMATIC_VARIANTS } from '../../modules/local/isec_vcfs'
//include { PCGR_VCF as PCGR_READY_VCF                      } from '../../modules/local/pcgr_vcf'

//This will run the PCGR
//include { PCGR as RUN_PCGR } from '../modules/local/pcgr'

//We'll leave CPSR for the future
//include { CPSR as RUN_CPSR       } from '../modules/local/cpsr'
//include { CPSR_VALIDATE_INPUT    } from '../../modules/local/validate_input'


//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nfpgcrdev_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFPGCRDEV {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    /*FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())*/

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'nfpgcrdev_software_'  + 'mqc_'  + 'versions.yml',
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
