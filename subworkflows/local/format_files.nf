/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF                        } from '../../modules/local/reformat_input/reformat_vcf'
include { REFORMAT_CNA                        } from '../../modules/local/reformat_input/reformat_cna'
include { INTERSECT_SOMATIC_VARIANTS          } from '../../modules/local/reformat_input/isec_vcf'
include { BCFTOOLS_ISEC                       } from '../../modules/nf-core/bcftools/isec/main'
include { PCGR_VCF                            } from '../../modules/local/reformat_input/pcgr_vcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

pcgr_header = Channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists:true)


workflow FORMAT_FILES {

    take:
    vcf_files
    cna_files

    main:
    ch_versions = Channel.empty()


    // Reformat input files
    REFORMAT_VCF( vcf_files )
    REFORMAT_CNA( cna_files )

    vcf_ch = REFORMAT_VCF.out.vcf
    cna_ch = REFORMAT_CNA.out.cna

    // vcf_ch.view()
    // cna_ch.view()


    // Intersect somatic variants
    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic      = vcf_ch.map{ meta, vcf, tbi -> var = [:]; var.patient = meta.patient; var.status = meta.status; var.sample = meta.sample; return [ var, vcf, tbi ] }
    per_sample_somatic_vcfs = per_sample_somatic.map{ meta, vcf, tbi -> return [ var, vcf, tbi ] }.groupTuple()

    /*per_sample_somatic_vcfs.map{ it ->
            "My values are:\n$it\n"
        }
        .collectFile(
            name: 'INTERSECT_SOMATIC_VARIANTS.txt',
            storeDir:'.',
            keepHeader: true,
            skip: 1
        )*/



    per_sample_somatic_vcfs = per_sample_somatic_vcfs.map{ 
        meta, vcfs, tbis ->
        [meta, vcfs, tbis, (1..vcfs.size()).toList()]
    }.branch {
        meta, vcfs, tbis, vcf_size ->
            single: vcf_size < 2
            multiple: vcf_size > 1
    }
    per_sample_somatic_vcfs_multiple = per_sample_somatic_vcfs.multiple.transpose(by: 3).map{
        meta, vcfs, tbis, vcf_size ->
        [ meta + ['vcf_size': vcf_size], vcfs, tbis ]
    }


    per_sample_somatic_vcfs_multiple.dump(tag: 'per_sample_somatic_vcfs')
    BCFTOOLS_ISEC ( per_sample_somatic_vcfs_multiple )
    
    // TODO : per_sample_somatic_vcfs.single have to go to BCFTOOLS VIEW process to retrieve create the sites.txt files
    
    BCFTOOLS_ISEC.out.results.dump(tag: 'BCFTOOLS_ISEC_out_results')


    ch_isec_somatic_postprocess = BCFTOOLS_ISEC.out.results.map {
        meta, results ->
        [ meta.subMap([ 'patient', 'status', 'sample']), results ]
    }.groupTuple()

    ch_isec_somatic_postprocess.dump(tag: 'ch_isec_somatic_postprocess')
    
    INTERSECT_SOMATIC_VARIANTS( ch_isec_somatic_postprocess ) 

    BCFTOOLS_ISEC ( per_sample_somatic_vcfs )

    INTERSECT_SOMATIC_VARIANTS( per_sample_somatic_vcfs )

    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs)

    /*sample_vcfs_keys.map{ it ->
            "My values are:\n$it\n"
        }
        .collectFile(
            name: 'PCGR_VCF.txt',
            storeDir:'.',
            keepHeader: true,
            skip: 1
        )*/

    PCGR_VCF( sample_vcfs_keys, pcgr_header.collect() )



    ch_versions = ch_versions.mix( REFORMAT_VCF.out.versions )
    ch_versions = ch_versions.mix( REFORMAT_CNA.out.versions )
    ch_versions = ch_versions.mix( INTERSECT_SOMATIC_VARIANTS.out.versions )
    ch_versions = ch_versions.mix( PCGR_VCF.out.versions )


    //PCGR_VCF.out.vcf.view()
    //cna_ch.view()

    emit:
    pcgr_ready_vcf = params.cna_analysis ? PCGR_VCF.out.vcf.join( cna_ch ) : PCGR_VCF.out.vcf.map{ meta, vcf, tbi -> return [ meta, vcf, tbi, [] ] }
    versions       = ch_versions               // channel: [ path(versions.yml) ]

}
