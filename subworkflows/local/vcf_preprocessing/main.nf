/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'

include { BCFTOOLS_NORM    } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_FILTER  } from '../../../modules/nf-core/bcftools/filter/main'
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
    ch_versions = channel.empty()


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

    TABIX_BGZIPTABIX.out.gz_tbi.set { ch_vcf_bgzipped }

    // Create index for files that need tabix
    ch_vcf.to_tabix
        .map { meta, vcf, _tbi ->
            vcf: [meta, vcf]
        }
        .set { ch_vcf_to_tabix }

    TABIX_TABIX(ch_vcf_to_tabix)

    ch_vcf_to_tabix
        .join(TABIX_TABIX.out.tbi)
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

    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

    emit:
    filtered_ch // channel: [ meta, path(vcf_file), path(tbi_file), path(cna_file) ]
    ch_cna_files
    versions     = ch_versions // channel: [ path(versions.yml) ]
}
