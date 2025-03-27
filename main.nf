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

include { VCF_PREPROCESSING       } from './workflows/vcf_preprocessing'
include { PCGR_PREPROCESSING      } from './workflows/pcgr_preprocessing'
include { PCGR_PIPELINE           } from './workflows/pgcr_pipeline'
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
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()


pcgr_header = Channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists:true)

if (params.database) { 
    pcgr_dir = Channel.fromPath("${params.database}/data/${params.genome.toLowerCase()}") 
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
workflow NF_VARIANT_PRIORITIZATION {

    take:
    samplesheet

    main:
    //samplesheet.view()
    //from samplesheet divide vcf from cna
    samplesheet.multiMap{ meta, vcf, tbi, cna ->
                                vcf_files : [ meta, vcf, tbi ]
                                cna_files : [ meta, cna ]
                                }.set{ files }

    // simplify metadata in CNA channels for re-merging to processed VCF files.
    vcf_files = files.vcf_files

    //vcf_files.view()

    ch_cna_files = files.cna_files.map{ 
        meta, cna -> var = [:]; var.patient = meta.patient; var.status = meta.status; var.sample = meta.sample;
        return [ meta, cna ] }
   
    //ch_cna_files.view()
   
    //
    // WORKFLOW: Run pipeline
    //
    VCF_PREPROCESSING(
            vcf_files,
            fasta
    )

    normalized_somatic = VCF_PREPROCESSING.out.normalised_somatic

    //normalized_somatic.view()

    PCGR_PREPROCESSING (
        normalized_somatic,
        ch_cna_files,
        pon_vcf,
        pcgr_header
    )

    preprocessed_vcf = PCGR_PREPROCESSING.out.pcgr_ready_vcf

    PCGR_PIPELINE (
        preprocessed_vcf,
        pon_vcf,
        vep_cache,
        pcgr_dir
    ) 

    //
    // WORKFLOW: Run pipeline
    //
    /*CPGR_PREPROCESSING (
        samplesheet
    ) 

    CPGR_PIPELINE (
        samplesheet
    )*/ 


    /*emit:
    multiqc_report = NFPGCRDEV.out.multiqc_report // channel: /path/to/multiqc_report.html*/
}


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

    //samplesheet = PIPELINE_INITIALISATION.out.samplesheet

    //
    // WORKFLOW: Run main workflow
    //

    NF_VARIANT_PRIORITIZATION (
        PIPELINE_INITIALISATION.out.samplesheet //,
        //fasta
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    /*PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NF_PCGR.out.multiqc_report
    )*/
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
