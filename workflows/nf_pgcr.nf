/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//This optionally tabix and bgzip the input VCF file [input_check]
//Questi due processi andranno valutate in pipeline initialization
//include { TABIX_TABIX as TABIX_INPUT_VCF      } from '../../modules/nf-core/tabix/tabix/main'
//include { TABIX_BGZIPTABIX as BGZIP_INPUT_VCF } from '../../modules/nf-core/tabix/bgziptabix/main'

//This reformat the input VCF file [format_files]
include { BCFTOOLS_NORM                    } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER                  } from '../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX as TABIX_FILTERED    } from '../../modules/nf-core/tabix/tabix/main'
include { REFORMAT_VCF                     } from '../../modules/local/pcgr_reformat'
include { REFORMAT_CNA                     } from '../../modules/local/pcgr_reformat'
include { REFORMAT_PON                     } from '../../modules/local/pcgr_reformat' 

//This merge the files [merge_vcfs]
include { BCFTOOLS_CONCAT                                 } from '../../modules/nf-core/bcftools/concat/main'
include { TABIX_TABIX as TABIX_CONCAT                     } from '../../modules/nf-core/tabix/tabix/main'
include { ISEC_SOMATIC_VCFS as INTERSECT_SOMATIC_VARIANTS } from '../../modules/local/isec_vcfs'
include { PCGR_VCF as PCGR_READY_VCF                      } from '../../modules/local/pcgr_vcf'

//This will run the PCGR
include { PCGR as RUN_PCGR                           } from '../modules/local/pcgr'

//We'll leave CPSR for the future
include { CPSR_VALIDATE_INPUT                      } from '../../modules/local/validate_input'
include { TABIX_BGZIPTABIX as BGZIPTABIX_CPSR      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { CPSR as RUN_CPSR                         } from '../modules/local/cpsr'

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
    fasta
    vcf_files
   
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    BCFTOOLS_NORM( vcf_files, fasta )
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    BCFTOOLS_FILTER( BCFTOOLS_NORM.out.vcf )
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())
    
    TABIX_FILTERED( BCFTOOLS_FILTER.out.vcf )
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(TABIX_FILTERED.out.versions.first())

    normalised_somatic  = BCFTOOLS_FILTER.out.vcf.join( TABIX_FILTERED.out.tbi ).filter{ it -> meta = it[0]; meta.status == 'somatic' }

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


workflow PCGR_PREPROCESSING {

    take:
    normalised_somatic
    cna_files
    pon_vcf
    pcgr_header
   
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //PCGR 
    REFORMAT_VCF( normalised_somatic )
    ch_versions = ch_versions.mix(REFORMAT_VCF.out.versions.first())

    REFORMAT_CNA( cna_files )
    ch_versions = ch_versions.mix(REFORMAT_CNA.out.versions.first())

    REFORMAT_PON( pon_vcf )
    ch_versions = ch_versions.mix(REFORMAT_PON.out.versions.first())

    copy_number = REFORMAT_CNA.out.cna.map{ meta, tsv -> var = [:]; var.id = meta.id; var.patient = meta.patient; var.status = meta.status; var.sample = meta.sample; var.tool = meta.tool; return [ var, tsv ] }

    somatic_files = params.cna_analysis ? REFORMAT_VCF.out.vcf.join( copy_number ) : REFORMAT_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

    pon = params.tumor_only && params.pon_vcf ? REFORMAT_PON.out.pon : pon_vcf

    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic = somatic_files.map{ meta, vcf, tbi, cna -> var = [:]; var.patient = meta.patient; var.sample = meta.sample; return [ var, vcf, tbi, cna ] }
    per_sample_somatic_vcfs = per_sample_somatic.map{ meta, vcf, tbi, cna -> return [ var, vcf, tbi ] }.groupTuple()
    INTERSECT_SOMATIC_VARIANTS( per_sample_somatic_vcfs )

    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs)
    PCGR_READY_VCF( sample_vcfs_keys, pcgr_header.collect() )
    pcgr_ready_vcf = params.cna_analysis ? PCGR_READY_VCF.out.vcf.join( per_sample_somatic.map{ meta, vcf, tbi, cna -> return [ meta, cna ] } ) : PCGR_READY_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

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

    emit:
    pcgr_ready_vcf 
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}



workflow PCGR_PIPELINE {

    take:
    pcgr_ready_vcf,
    pon_vcf,
    vep_cache,
    pcgr_dir
   
    main:

    ch_versions = Channel.empty()
 

    //PCGR
    RUN_PCGR( pcgr_ready_vcf, ch_pcgr_dir.collect(), pon_vcf, vep_cache.collect()  )

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



}



/*workflow CPGR_PREPROCESSING {

    take:
    fasta
    vcf_files
    cna_files
    pon_vcf
    pcgr_header
   
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //CPGR
    normalised_germline = BCFTOOLS_FILTER.out.vcf.join( TABIX_FILTERED.out.tbi ).filter{ it -> meta = it[0]; meta.status == 'germline' }
    //Aggiungere i moduli di CPGR
    CPSR_VALIDATE_INPUT( germline_files, pcgr_dir.collect() )
    BGZIPTABIX_CPSR( CPSR_VALIDATE_INPUT.out.validated_vcf )
    per_sample_germline = BGZIPTABIX_CPSR.out.gz_tbi.map{ meta, vcf, tbi -> var = [:]; var.patient = meta.patient; var.sample = meta.sample; var.id = "${meta.patient}.${meta.sample}"; return [ var, vcf, tbi ] }.groupTuple()
    BCFTOOLS_CONCAT( per_sample_germline )
    TABIX_CONCAT( BCFTOOLS_CONCAT.out.vcf )
    cpsr_ready_vcf = BCFTOOLS_CONCAT.out.vcf.join( TABIX_CONCAT.out.tbi )

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

}



workflow PCGR_PIPELINE {

    take:
    pcgr_ready_vcf
    pcgr_dir
   
    main:

    ch_versions = Channel.empty()
 

    //PCGR
    RUN_PCGR( pcgr_ready_vcf, ch_pcgr_dir.collect(), pon_vcf, vep_cache.collect()  )

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



}*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
