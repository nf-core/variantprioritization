/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_CONCAT   } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_REHEADER } from '../../../modules/nf-core/bcftools/reheader'
include { CREATESAMPLEFILE  } from '../../../modules/local/createsamplefile'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_GERMLINE {
    take:
    normalised_germline

    main:
    // Reheader and combine germline VCFs
    CREATESAMPLEFILE(normalised_germline.map { meta, _vcf, _tbi -> meta })

    ch_vcf_sample = normalised_germline
        .join(CREATESAMPLEFILE.out.samplefile)
        .map { meta, vcf, _tbi, samplefile -> [meta, vcf, [], samplefile] }

    BCFTOOLS_REHEADER(
        ch_vcf_sample,
        [[], []],
    )

    ch_vcf_index_rh = BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index)

    ch_concat_in = ch_vcf_index_rh
        .map { meta, vcf, tbi -> tuple(meta.sample, meta, vcf, tbi) }
        .groupTuple(by: 0)
        .map { _sample, metas, vcfs, tbis ->
            def id = metas[0].sample
            def sample = metas[0].sample
            def patient = metas[0].patient
            def status = metas[0].status
            def tools = metas.collect{ maps -> maps.tool }.join(",")
            def new_meta = [
                id: id,
                patient: patient,
                sample: sample,
                status: status,
                tool: tools
            ]
            tuple(new_meta, vcfs, tbis)
        }

    BCFTOOLS_CONCAT(ch_concat_in)

    combined_germline = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi)

    emit:
    combined_germline
}
