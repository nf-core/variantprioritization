/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_CONCAT   } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_FILTER   } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_NORM     } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER } from '../../../modules/nf-core/bcftools/reheader'
include { CREATESAMPLEFILE  } from '../../../modules/local/createsamplefile'
include { TABIX_BGZIPTABIX  } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_PREPROCESSING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input and processed by PIPELINE_INITIALISATION
    fasta

    main:
    ch_samplesheet
        .multiMap { meta, vcf, tbi, cna ->
            vcf_files: [meta, vcf, tbi]
            cna_files: [meta, cna]
        }
        .set { files }

    // simplify metadata in CNA channels for re-merging to processed VCF files.
    vcf_files = files.vcf_files
    ch_cna_files = files.cna_files
        .map { meta, cna ->
            def var = [:]
            var.patient = meta.patient
            var.status = meta.status
            var.sample = meta.sample
            return [var, cna]
        }
        .distinct()

    // Create subchannels for files that need bgzipping and tabix indexing
    ch_vcf = vcf_files.branch { items ->
        to_bgzip: items[0].bgzip_vcf == false
        to_tabix: items[0].tabix_vcf == false
        ready: true
    }

    // Process files that need bgzipping
    ch_vcf.to_bgzip
        .map { meta, vcf, _tbi ->
            vcf: [meta, vcf]
        }
        .set { ch_vcf_to_bgzip }

    TABIX_BGZIPTABIX(ch_vcf_to_bgzip)

    TABIX_BGZIPTABIX.out.gz_index.set { ch_vcf_bgzipped }

    // Create index for files that need tabix
    ch_vcf.to_tabix
        .map { meta, vcf, _tbi -> tuple(meta.id, meta, vcf) }
        .set { ch_vcf_to_tabix_keyed }

    ch_vcf_to_tabix_keyed
        .map { _id, meta, vcf -> [meta, vcf] }
        .set { ch_vcf_to_tabix }

    TABIX_TABIX(ch_vcf_to_tabix)

    ch_vcf_to_tabix_keyed
        .join(
            TABIX_TABIX.out.index.map { meta, tbi -> tuple(meta.id, meta, tbi) },
            by: 0
        )
        .map { _id, meta, vcf, _meta2, tbi -> [meta, vcf, tbi] }
        .set { ch_vcf_with_tabix }

    // Combine all processed files into a final channel
    ch_vcf.ready
        .mix(ch_vcf_bgzipped)
        .mix(ch_vcf_with_tabix)
        .set { ch_vcf }

    //VCF_PREPROCESSING
    BCFTOOLS_NORM(ch_vcf, fasta)
    norm_ch = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi)

    BCFTOOLS_FILTER(norm_ch)
    filtered_ch = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi)

    normalised_germline = filtered_ch.filter { meta, _vcf, _tbi -> meta.status == 'germline' }
    normalised_somatic = filtered_ch.filter { meta, _vcf, _tbi -> meta.status == 'somatic' }


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

    ch_vcf_index_rh.dump(tag: 'ch_vcf_index_rh')

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

    combined_germline.dump(tag: 'combined_germline')

    emit:
    combined_germline
    normalised_somatic
    ch_cna_files
}
