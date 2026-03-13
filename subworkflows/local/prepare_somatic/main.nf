/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF    } from '../../../modules/local/reformat/vcf'
include { REFORMAT_CNA    } from '../../../modules/local/reformat/cna'
include { INTERSECT_VCF   } from '../../../modules/local/intersect/vcf'
include { PCGR_PREPAREVCF } from '../../../modules/local/pcgr/preparevcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_SOMATIC {
    take:
    vcf_files
    cna_files

    main:
    pcgr_header = channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists: true)


    // Reformat input files
    REFORMAT_VCF(vcf_files)
    REFORMAT_CNA(cna_files)

    vcf_ch = REFORMAT_VCF.out.vcf
    cna_ch = REFORMAT_CNA.out.cna

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

    per_sample_somatic_vcfs = per_sample_somatic
        .map { var, vcf, tbi ->
            return [var, vcf, tbi]
        }
        .groupTuple()

    INTERSECT_VCF(per_sample_somatic_vcfs)

    sample_vcfs_keys = INTERSECT_VCF.out.variant_tool_map
        .join(per_sample_somatic_vcfs)
        .map { meta, keys, _meta2, vcfs, tbis -> [meta, keys, vcfs, tbis] }

    PCGR_PREPAREVCF(sample_vcfs_keys, pcgr_header.collect())

    emit:
    pcgr_ready_vcf = params.cna_analysis
        ? PCGR_PREPAREVCF.out.vcf.join(cna_ch)
        : PCGR_PREPAREVCF.out.vcf.map { meta, vcf, tbi ->
            return [meta, vcf, tbi, []]
        }
}
