/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_bwa_align_stats_pipeline'
include { FASTP                     } from '../modules/nf-core/fastp'
include { FASTA_INDEX_DNA           } from '../subworkflows/nf-core/fasta_index_dna'
include { FASTQ_ALIGN_DNA           } from '../subworkflows/nf-core/fastq_align_dna'
include { BAM_SORT_STATS_SAMTOOLS   } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_MPILEUP_SNPS_VCF      } from '../subworkflows/local/bam_mpileup_snps_vcf'  

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

    // 
    // QC of input reads
    //

    FASTP(
        ch_samplesheet,
        [], // we are not using any adapter fastas at the moment
        false, // we don't use discard_trimmed_pass
        params.fastp_save_trimmed_fail,
        params.fastp_save_merged
    )

    // 
    // Create index for reference genome, empty liftover channel
    //

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

    FASTQ_ALIGN_DNA(
        FASTP.out.reads,
        FASTA_INDEX_DNA.out.index.first(),
        ch_refgenome,
        params.aligner,
        params.sort_bam )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

    // Sort resulting bam and produce al

    BAM_SORT_STATS_SAMTOOLS(
        FASTQ_ALIGN_DNA.out.bam,
        ch_refgenome    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // optional: produce pileup file

    if ( !params.skip_pileup) {
        BAM_MPILEUP_SNPS_VCF( 
            BAM_SORT_STATS_SAMTOOLS.out.bam,
            ch_refgenome)
        ch_versions = ch_versions.mix(BAM_MPILEUP_SNPS_VCF.out.versions)
    }

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
    index          = FASTA_INDEX_DNA.out.index              // channel: [ val(meta), path(bwamem2_dir)]
    versions       = ch_versions                            // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
