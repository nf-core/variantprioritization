/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//This reformat the input VCF file [format_files]
include { REFORMAT_VCF                     } from '../modules/local/pcgr_reformat'
include { REFORMAT_CNA                     } from '../modules/local/pcgr_reformat'
include { REFORMAT_PON                     } from '../modules/local/pcgr_reformat' 

//This merge the files [merge_vcfs]
include { BCFTOOLS_CONCAT                                 } from '../modules/nf-core/bcftools/concat/main'
include { TABIX_TABIX as TABIX_CONCAT                     } from '../modules/nf-core/tabix/tabix/main'
include { ISEC_SOMATIC_VCFS as INTERSECT_SOMATIC_VARIANTS } from '../modules/local/isec_vcfs'
include { PCGR_VCF as PCGR_READY_VCF                      } from '../modules/local/pcgr_vcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    //REFORMAT_VCF.out.vcf.view()
    //REFORMAT_CNA.out.cna.view()

    /*copy_number = REFORMAT_CNA.out.cna.map{ meta, tsv ->    var = [:];
                                                            var.id = meta.id;
                                                            var.patient = meta.patient;
                                                            var.status = meta.status;
                                                            var.sample = meta.sample;
                                                            var.tool = meta.tool;
                                                            return [ var, tsv ] }*/

    //copy_number.view()

    somatic_files = params.cna_analysis ? REFORMAT_VCF.out.vcf.join( REFORMAT_CNA.out.cna ) 
                                        : REFORMAT_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

    somatic_files.view()
    
    pon = params.tumor_only && params.pon_vcf ? REFORMAT_PON.out.pon : pon_vcf

    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]

    per_sample_somatic = somatic_files.map{ meta, vcf, tbi, cna ->    var = [:];
                                                                        var.patient = meta.patient;
                                                                        var.sample = meta.sample;
                                                                        return [ var, vcf, tbi, cna ] }

    per_sample_somatic_vcfs = per_sample_somatic.map{ meta, vcf, tbi, cna -> return [ var, vcf, tbi ] }
                                                .groupTuple()


    INTERSECT_SOMATIC_VARIANTS( per_sample_somatic_vcfs )

    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs)
    PCGR_READY_VCF( sample_vcfs_keys, pcgr_header.collect() )
    pcgr_ready_vcf = params.cna_analysis ? PCGR_READY_VCF.out.vcf.join( per_sample_somatic.map{ meta, vcf, tbi, cna -> return [ meta, cna ] } ) : PCGR_READY_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }

    //
    // Collate and save software versions
    //
    /*softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'nfpgcrdev_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }*/

    emit:
    pcgr_ready_vcf 
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
