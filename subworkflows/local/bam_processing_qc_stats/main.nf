
include { SAMTOOLS_INDEX                    } from '../../../modules/nf-core/samtools/index/main'
include { GATK4SPARK_MARKDUPLICATES         } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { SAMTOOLS_STATS                    } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH                          } from '../../../modules/nf-core/mosdepth/main'

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
    ch_reports  = Channel.empty()

    // Mark duplicates and dort resulting bam file
    GATK4SPARK_MARKDUPLICATES(
        ch_bam,
        ch_ref_fasta.map{ meta, fasta -> fasta },
        ch_fasta_fai.map{ meta, fasta_fai -> fasta_fai },
        ch_dictionary.map{ meta, dict -> dict}
    )
    SAMTOOLS_INDEX( GATK4SPARK_MARKDUPLICATES.out.output )
    
    // Join bam + bai, calculate mapping stats and depth
    GATK4SPARK_MARKDUPLICATES.out.output
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch: true )
        .set{ ch_bam_bai }
    SAMTOOLS_STATS(
        ch_bam_bai, 
        ch_ref_fasta
        )
    if (!params.skip_coverage){
        MOSDEPTH(
            ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, bai, [] ]},
            ch_ref_fasta
        )
    }

    // Gather all reports generated
    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)
    ch_reports = ch_reports.mix(MOSDEPTH.out.global_txt)
    ch_reports = ch_reports.mix(MOSDEPTH.out.regions_txt)

    // Capture versions
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_versions = ch_versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    bam      = GATK4SPARK_MARKDUPLICATES.out.output     // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai                   // channel: [ val(meta), [ bai ] ]
    bam_bai  = ch_bam_bai                               // channel: [ val(meta), [ bam ], [bai] ]
    csi      = SAMTOOLS_INDEX.out.csi                   // channel: [ val(meta), [ csi ] ]
    reports  = ch_reports                               // channel: [ val(meta), [report1.txt], [report2.txt], ... ]   
    versions = ch_versions                              // channel: [ versions.yml ]
}

