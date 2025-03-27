/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//This optionally tabix and bgzip the input VCF file [input_check]
include { TABIX_TABIX as TABIX_INPUT_VCF      } from '../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as BGZIP_INPUT_VCF } from '../modules/nf-core/tabix/bgziptabix/main'

//This reformat the input VCF file [format_files]
include { BCFTOOLS_NORM                    } from '../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_NORM        } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER                  } from '../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX as TABIX_FILTERED    } from '../modules/nf-core/tabix/tabix/main'

//multiqc
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

workflow VCF_PREPROCESSING {

    take:
    vcf_files
    fasta
   
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //vcf_files.view()

    // Assuming your input channel is named 'input_channel'
    process_channel = vcf_files.branch {
        to_bgzip: it[0].bgzip_vcf == true
        to_tabix: it[0].tabix_vcf == true
        already_processed: it[0].bgzip_vcf == false && it[0].tabix_vcf == false
    }

    // Process files that need bgzipping
    /*process_channel.to_bgzip.view()
    BGZIP_INPUT_VCF(process_channel.to_bgzip)
    bgzipped_channel = BGZIP_INPUT_VCF.out.gz_tbi
    //bgzipped_channel.view()*/

    /*
    // Create index for files that need tabix
    //process_channel.to_tabix.view()
    TABIX_INPUT_VCF(process_channel.to_tabix)
    tabix_channel = TABIX_INPUT_VCF.out.tbi */

    /*tabixed_channel = process_channel.to_tabix.join(tabix_channel)


    //tabixed_channel.view()*/

    // Combine all processed files into a final channel
    //vcf_files.view()

    //fasta.view()
    BCFTOOLS_NORM( 
        vcf_files,
        fasta 
        )

    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())


    TABIX_NORM( BCFTOOLS_NORM.out.vcf)

    //BCFTOOLS_NORM.out.tbi.view()
    tabixed_channel = BCFTOOLS_NORM.out.vcf.join(TABIX_NORM.out.tbi)
    //tabixed_channel.view()

    BCFTOOLS_FILTER( tabixed_channel )
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    
    TABIX_FILTERED( BCFTOOLS_FILTER.out.vcf )
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(TABIX_FILTERED.out.versions.first())

    normalised_somatic  = BCFTOOLS_FILTER.out.vcf.join( TABIX_FILTERED.out.tbi ).filter{ it -> meta = it[0]; meta.status == 'somatic' }
   
    //normalised_somatic.view()
   
    //normalised_germline  = BCFTOOLS_FILTER.out.vcf.join( TABIX_FILTERED.out.tbi ).filter{ it -> meta = it[0]; meta.status == 'somatic' }
    //normalised_germline.view()

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
    /*ch_multiqc_config        = Channel.fromPath(
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
    )*/

    emit:
    //multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    normalised_somatic
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
