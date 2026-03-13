/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER  } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_NORM    } from '../../../modules/nf-core/bcftools/norm'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_PREPROCESSING {
    take:
    ch_samplesheet
    fasta

    main:
    ch_vcf_input = ch_samplesheet.map { meta, vcf, tbi, _cna -> [meta, vcf, tbi] }

    ch_cna_files = ch_samplesheet
        .map { meta, _vcf, _tbi, cna ->
            [meta.subMap(['patient', 'status', 'sample']), cna]
        }
        .distinct()

    ch_vcf_branched = ch_vcf_input.branch { items ->
        to_bgzip: items[0].bgzip_vcf == false
        to_tabix: items[0].tabix_vcf == false
        ready: true
    }

    ch_vcf_for_bgzip = ch_vcf_branched.to_bgzip.map { meta, vcf, _tbi -> [meta, vcf] }

    TABIX_BGZIPTABIX(ch_vcf_for_bgzip)

    TABIX_BGZIPTABIX.out.gz_index.set { ch_vcf_bgzipped }

    ch_vcf_for_tabix = ch_vcf_branched.to_tabix.map { meta, vcf, _tbi -> [meta, vcf] }
    ch_vcf_for_tabix_keyed = ch_vcf_for_tabix.map { meta, vcf -> tuple(meta.id, meta, vcf) }

    TABIX_TABIX(ch_vcf_for_tabix)

    ch_vcf_for_tabix_keyed
        .join(
            TABIX_TABIX.out.index.map { meta, tbi -> tuple(meta.id, tbi) },
            by: 0
        )
        .map { _id, meta, vcf, tbi -> [meta, vcf, tbi] }
        .set { ch_vcf_with_tabix }

    ch_vcf_branched.ready
        .mix(ch_vcf_bgzipped)
        .mix(ch_vcf_with_tabix)
        .set { ch_vcf_preprocessed }

    BCFTOOLS_NORM(ch_vcf_preprocessed, fasta)
    ch_norm = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi)

    BCFTOOLS_FILTER(ch_norm)
    ch_filtered = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)

    normalised_germline = ch_filtered.filter { meta, _vcf, _tbi -> meta.status == 'germline' }
    normalised_somatic = ch_filtered.filter { meta, _vcf, _tbi -> meta.status == 'somatic' }

    emit:
    normalised_germline
    normalised_somatic
    ch_cna_files
}
