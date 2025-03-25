#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/nfpgcrdev
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/nfpgcrdev
    Website: https://nf-co.re/nfpgcrdev
    Slack  : https://nfcore.slack.com/channels/nfpgcrdev
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VCF_PREPROCESSING       } from './workflows/nf_pgcr'
include { PCGR_PREPROCESSING      } from './workflows/nf_pgcr'
include { PCGR_PIPELINE           } from './workflows/nf_pgcr'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_nfpgcrdev_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_nfpgcrdev_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_nfpgcrdev_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`

params.fasta = getGenomeAttribute('fasta')

ch_fasta    = Channel.fromPath(params.fasta, checkIfExists: true)

pcgr_header = Channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists:true)

if (params.database) { 
    ch_pcgr_dir = Channel.fromPath("${params.database}/data/${params.genome.toLowerCase()}") 
    } 
    else 
    { exit 1, "Please provide a path to the PCGR annotation database." }

if (params.tumor_only && params.pon_vcf) { 
    pon_vcf = file(params.pon_vcf, checkIfExists: true) 
    } else 
    { pon_vcf = Channel.value([]) }

vep_cache               = Channel.fromPath(params.vep_cache        )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NF_PCGR {

    take:
    samplesheet, // channel: samplesheet read in from --input

    main:

    //from samplesheet divide vcf from cna

    //
    // WORKFLOW: Run pipeline
    //
    VCF_PREPROCESSING(
            fasta,
            vcf_files
    )

    normalized_vcf = VCF_PREPROCESSING.out.normalised_somatic

    PGCR_PREPROCESSING (
        normalized_vcf,
        cna_files,
        pon_vcf,
        pcgr_header
    )

    preprocessed_vcf = PGCR_PREPROCESSING.out.pcgr_ready_vcf

    PCGR_PIPELINE (
        preprocessed_vcf,
        pon_vcf,
        vep_cache,
        pcgr_dir
    )  

    emit:
    multiqc_report = NFPGCRDEV.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*workflow NF_CPGR {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CPGR_PREPROCESSING (
        samplesheet
    ) 

    NF_CPGR (
        samplesheet
    ) 


    emit:
    multiqc_report = NFPGCRDEV.out.multiqc_report // channel: /path/to/multiqc_report.html
}*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input //samplesheet
    )

    //
    // WORKFLOW: Run main workflow
    //
    NF_PCGR (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NF_NFPGCRDEV.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
