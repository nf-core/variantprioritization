/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF               } from '../../../modules/local/reformat/vcf'
include { REFORMAT_CNA               } from '../../../modules/local/reformat/cna'
include { INTERSECT_SOMATIC_VARIANTS } from '../../../modules/local/isec_vcf'
include { PCGR_PREPAREVCF            } from '../../../modules/local/pcgr/preparevcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FORMAT_FILES {
    take:
    vcf_files
    cna_files

    main:
    ch_versions = channel.empty()

    pcgr_header = channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists: true)


    // Reformat input files
    REFORMAT_VCF(vcf_files)
    REFORMAT_CNA(cna_files)

    vcf_ch = REFORMAT_VCF.out.vcf
    cna_ch = REFORMAT_CNA.out.cna

    // vcf_ch.view()
    // cna_ch.view()


    // Intersect somatic variants
    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic = vcf_ch.map { meta, vcf, tbi ->
        def var = [:]
        var.patient = meta.patient
        var.status = meta.status
        var.sample = meta.sample
        return [var, vcf, tbi]
    }
    per_sample_somatic_vcfs = per_sample_somatic.map { var, vcf, tbi ->
        return [var, vcf, tbi]
    }.groupTuple()

    /*per_sample_somatic_vcfs.map{ it ->
            "My values are:\n$it\n"
        }
        .collectFile(
            name: 'INTERSECT_SOMATIC_VARIANTS.txt',
            storeDir:'.',
            keepHeader: true,
            skip: 1
        )*/

    INTERSECT_SOMATIC_VARIANTS(per_sample_somatic_vcfs)

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

    PCGR_PREPAREVCF(sample_vcfs_keys, pcgr_header.collect())

    ch_versions = ch_versions.mix(INTERSECT_SOMATIC_VARIANTS.out.versions)
    ch_versions = ch_versions.mix(PCGR_PREPAREVCF.out.versions)

    emit:
    pcgr_ready_vcf = params.cna_analysis ? PCGR_PREPAREVCF.out.vcf.join(cna_ch) : PCGR_PREPAREVCF.out.vcf.map { meta, vcf, tbi ->
        return [meta, vcf, tbi, []]
    }
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
