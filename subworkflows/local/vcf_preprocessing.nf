/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX      } from '../../modules/nf-core/tabix/tabix/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_PREPROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input and processed by PIPELINE_INITIALISATION
   
    main:
    ch_versions = Channel.empty()

    // Create subchannels for files that need bgzipping and tabix indexing

    ch_vcf = ch_samplesheet.branch {
        to_bgzip: it[0].bgzip_vcf == false
        to_tabix: it[0].tabix_vcf == false
        ready: true
    }

    // Process files that need bgzipping

    ch_vcf.to_bgzip
        .multiMap { meta, vcf, _tbi, cna ->
            vcf: [meta, vcf]
            cna: [meta, cna]
        }
        .set { ch_vcf_to_bgzip }

    TABIX_BGZIPTABIX( ch_vcf_to_bgzip.vcf )

    TABIX_BGZIPTABIX.out.gz_tbi
        .join(ch_vcf_to_bgzip.cna)
        .set{ ch_vcf_bgzipped }

    // Create index for files that need tabix

    ch_vcf.to_tabix
        .multiMap { meta, vcf, _tbi, cna ->
            vcf: [meta, vcf]
            cna: [meta, cna]
        }
        .set { ch_vcf_to_tabix }

    TABIX_TABIX( ch_vcf_to_tabix.vcf )

    ch_vcf_to_tabix.vcf
        .join(TABIX_TABIX.out.tbi)
        .join(ch_vcf_to_tabix.cna)
        .set{ ch_vcf_with_tabix }

    // Combine all processed files into a final channel

    ch_vcf.ready
        .mix(ch_vcf_bgzipped)
        .mix(ch_vcf_with_tabix)
        .set { ch_vcf }

    emit:
    ch_vcf                                         // channel: [ meta, path(vcf_file), path(tbi_file), path(cna_file) ]
    versions       = ch_versions                   // channel: [ path(versions.yml) ]

}