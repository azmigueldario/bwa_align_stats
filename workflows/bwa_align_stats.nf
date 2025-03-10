/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_bwa_align_stats_pipeline'

include { PREPARE_INPUTS            } from '../subworkflows/local/prepare_inputs'
include { FASTQ_ALIGN_DNA           } from '../subworkflows/nf-core/fastq_align_dna'
include { BAM_SORT_STATS_SAMTOOLS   } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { BAM_PROCESSING_QC_STATS   } from '../subworkflows/local/bam_processing_qc_stats'  
include { BAM_MPILEUP_SNPS_VCF      } from '../subworkflows/local/bam_mpileup_snps_vcf'  

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BWA_ALIGN_STATS {

    take:
    ch_fastq_reads              // tuple [meta, [fastq1, fastq2] ] from '--input samplesheet.csv'
    ch_refgenome                // tuple [meta, fasta] from '--ref_genome fasta'
    ref_genome_path             // params.ref_genome
    fastp_save_trimmed_fail     // params.fastp_save_trimmed_fail
    fastp_save_merged           // params.fastp_save_merged
    aligner                     // params.aligner
    sort_bam                    // params.sort_bam

    main:

    ch_versions = Channel.empty()

    // 
    // FastQ processing, prepare inputs
    //

    PREPARE_INPUTS(
        ch_fastq_reads,
        ch_refgenome,
        ref_genome_path,
        fastp_save_trimmed_fail,
        fastp_save_merged
    )

    //
    // Align reads to reference genome
    //

    FASTQ_ALIGN_DNA(
        PREPARE_INPUTS.out.reads,
        PREPARE_INPUTS.out.bwamem2_index,
        PREPARE_INPUTS.out.ref_fasta,
        aligner,
        sort_bam )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

    // 
    // Processing alignment file
    //


    /* Sort resulting bam and produce alignment stats

    BAM_SORT_STATS_SAMTOOLS(
        FASTQ_ALIGN_DNA.out.bam,
        ch_refgenome    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    //
    // Optional: produce pileup and variant calling files
    //      Skip with '--skip_pileup', or produce only the 
    //      mpileup without variant calling using '--skip_variants'
    //

    if ( !params.skip_pileup) {
        BAM_MPILEUP_SNPS_VCF( 
            BAM_SORT_STATS_SAMTOOLS.out.bam,
            BAM_SORT_STATS_SAMTOOLS.out.bai,
            ch_refgenome,
            params.save_mpileup)
        ch_versions = ch_versions.mix(BAM_MPILEUP_SNPS_VCF.out.versions)
    }
    */

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
    index          = PREPARE_INPUTS.out.bwamem2_index       // channel: [ val(meta), path(bwamem2_dir)]
    versions       = ch_collated_versions                   // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
