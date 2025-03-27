/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//We'll leave CPSR for the future
include { CPSR_VALIDATE_INPUT                      } from '../../modules/local/validate_input'
include { TABIX_BGZIPTABIX as BGZIPTABIX_CPSR      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { CPSR as RUN_CPSR                         } from '../modules/local/cpsr'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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

}*/

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

}*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
