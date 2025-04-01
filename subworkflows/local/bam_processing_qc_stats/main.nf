
include { SAMTOOLS_INDEX                    } from '../../../modules/nf-core/samtools/index/main'
include { GATK4SPARK_MARKDUPLICATES         } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { SAMTOOLS_STATS                    } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH                          } from '../../../modules/nf-core/mosdepth/main'
include { SAMTOOLS_REINDEX_BAM              } from '../../../modules/local/samtools/reindexbam/main'
include { GOLEFT_INDEXCOV                   } from '../../../modules/nf-core/goleft/indexcov/main'

workflow BAM_PROCESSING_QC_STATS {

    take:
    ch_bam                          // channel: [ val(meta), [ bam ] ]
    ch_ref_fasta                    // channel: [ val(meta), [ genome.fa ] ]
    ch_fasta_fai                    // channel: [ val(meta), [ genome.fasta.fai ] ]
    ch_dictionary                   // channel: [ val(meta), [*.dict] ]

    main:

    ch_versions     = Channel.empty()
    ch_reports      = Channel.empty()
    indexcov_out    = Channel.empty()

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

    // Alignment stats after deduplication
    SAMTOOLS_STATS(
        ch_bam_bai, 
        ch_ref_fasta
        )
    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)
    
    if (!params.skip_coverage){
        MOSDEPTH(
            ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, bai, [] ]},
            ch_ref_fasta
        )

        // clean bam index and eval coverage
        SAMTOOLS_REINDEX_BAM(ch_bam_bai, ch_ref_fasta, ch_fasta_fai)
        SAMTOOLS_REINDEX_BAM.out.output
            .map{[[id:"indexcov"], it[1], it[2]]}
            .groupTuple()
            .set{indexcov_ch}
        GOLEFT_INDEXCOV(indexcov_ch, ch_fasta_fai)

        ch_reports      = ch_reports.mix(MOSDEPTH.out.global_txt)
        ch_reports      = ch_reports.mix(MOSDEPTH.out.regions_txt)
        indexcov_out    = indexcov_out.mix(GOLEFT_INDEXCOV.out.output)

        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    }

    // Capture versions
    ch_versions = ch_versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    bam         = GATK4SPARK_MARKDUPLICATES.out.output     // channel: [ val(meta), [ bam ] ]
    bai         = SAMTOOLS_INDEX.out.bai                   // channel: [ val(meta), [ bai ] ]
    bam_bai     = ch_bam_bai                               // channel: [ val(meta), [ bam ], [bai] ]
    csi         = SAMTOOLS_INDEX.out.csi                   // channel: [ val(meta), [ csi ] ]
    reports     = ch_reports                               // channel: [ val(meta), [report1.txt], [report2.txt], ... ]   
    versions    = ch_versions                              // channel: [ versions.yml ]
    indexcov    = indexcov_out                             // channel: [ val(meta), [results_dir] ]   
}

