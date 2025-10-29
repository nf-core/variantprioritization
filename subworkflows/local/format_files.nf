/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF                        } from '../../modules/local/reformat_input/reformat_vcf'
include { REFORMAT_CNA                        } from '../../modules/local/reformat_input/reformat_cna'
include { INTERSECT_SOMATIC_VARIANTS          } from '../../modules/local/reformat_input/isec_vcf'
include { PCGR_VCF                            } from '../../modules/local/reformat_input/pcgr_vcf'
include { BCFTOOLS_ISEC                       } from '../../modules/nf-core/bcftools/isec/main'
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

    // This could be refactored with the .collect() groovy list method.
    per_sample_somatic_vcfs.transpose().map {
        meta, vcf, tbis ->
        def tool_name = vcf.toString().tokenize('.')[2]
        [ meta , tool_name, vcf, tbis ]
    }.groupTuple().map {
        meta, tool_names, vcfs, tbis ->
        [ meta + ['tools': tool_names] , vcfs, tbis ]
    }.set { per_sample_somatic_vcfs } 
    
    per_sample_somatic_vcfs.map{ 
        meta, vcfs, tbis ->
        [meta, vcfs, tbis, (1..vcfs.size()).toList()]
    }.branch {
        meta, vcfs, tbis, vcf_size ->
            single: vcfs.size() < 2
            multiple: vcfs.size() > 1
    }.set { per_sample_somatic_vcfs }

    per_sample_somatic_vcfs.multiple.transpose(by: 3).map{
        meta, vcfs, tbis, vcf_size ->
        [ meta + ['vcf_size': vcf_size], vcfs, tbis ]
    }.set { per_sample_somatic_vcfs_multiple }


    BCFTOOLS_ISEC ( per_sample_somatic_vcfs_multiple )
    
    ch_isec_somatic_postprocess = BCFTOOLS_ISEC.out.results.map {
        meta, results ->
        [ meta.subMap([ 'patient', 'status', 'sample', 'tools']), results ]
    }.groupTuple()

    
    INTERSECT_SOMATIC_VARIANTS( ch_isec_somatic_postprocess ) 



    // merge mapping key back with sample VCFs, produce PCGR ready VCFs.
    sample_vcfs_keys = INTERSECT_SOMATIC_VARIANTS.out.variant_tool_map.join(per_sample_somatic_vcfs_multiple)

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
