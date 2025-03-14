
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'

include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_COVERAGE     } from '../../../modules/nf-core/samtools/coverage/main'  
include { SAMTOOLS_MPILEUP      } from '../../../modules/nf-core/samtools/mpileup/main'
include { BCFTOOLS_MPILEUP      } from '../../../modules/nf-core/bcftools/mpileup/main'      

workflow BAM_PROCESSING_QC_STATS {

    take:
    ch_bam                          // channel: [ val(meta), [ bam ] ]
    ch_ref_fasta                    // channel: [ val(meta), [ genome.fa ] ]
    ch_fasta_fai                    // channel: [ val(meta), [ genome.fasta.fai ] ]
    ch_dictionary                   // channel: [ val(meta), [*.dict] ]

    main:

    ch_versions = Channel.empty()

        // Sort the bam alignment
    SAMTOOLS_INDEX ( ch_bam )

    // Capture versions
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // Join bam and bai
    ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch: true)


    emit:
    // TODO nf-core: edit emitted channels
    bam      = ch_bam                           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai           // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi           // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

