/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bwa_align_stats_pipeline'
include { FASTA_INDEX_DNA        } from '../subworkflows/nf-core/fasta_index_dna'
include { FASTQ_ALIGN_DNA        } from '../subworkflows/nf-core/fastq_align_dna'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BWA_ALIGN_STATS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_refgenome   // required fasta input

    main:

    ch_versions = Channel.empty()

    // Create index for reference genome

    ch_refgenome
        .map {meta, fasta -> [meta, ""]}
        .set{ch_dummy_altliftover}

    FASTA_INDEX_DNA(
        ch_refgenome,
            // provide empty tuple of [meta, fasta] instead of liftover
        ch_dummy_altliftover,
        params.aligner )
    ch_versions = ch_versions.mix(FASTA_INDEX_DNA.out.versions)

    // Align reads to reference genome



    FASTA_INDEX_DNA.out.index.view()
    ch_refgenome.view()
    println "$params.aligner"

    FASTQ_ALIGN_DNA(
        ch_samplesheet,
        FASTA_INDEX_DNA.out.index.first(),
        ch_refgenome,
        params.aligner,
        params.sort_bam )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'bwa_align_stats_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
