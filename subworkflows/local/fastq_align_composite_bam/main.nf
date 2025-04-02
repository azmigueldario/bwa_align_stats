include { COMBINE_GENOMES                   } from '../../../modules/local/combinegenomes/main.nf'
include { BWAMEM2_INDEX as COMPOSITE_INDEX  } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as COMPOSITE_ALIGN    } from '../../../modules/nf-core/bwamem2/mem/main'
include { EXTRACT_MAPPING_STATS             } from '../../../modules/local/extractmappingstats/main.nf'

workflow FASTQ_ALIGN_COMPOSITE_BAM {

    take:
    ch_reads        // channel: [mandatory] meta, reads
    genome_list     // params.genome_list_composite

    main:

    ch_versions = Channel.empty()

    // Check input reference genomes and create input channel
    if (genome_list && genome_list instanceof List && genome_list.size() == 2 ) {

        Channel
            .fromList(genome_list)
            .collect()
            .map {
                [ [id:'input_genomes_composite'], [it[0], it[1]] ]
            }
            .set{input_composite_ch}
        input_composite_ch.view()
    }
    else { 
        error "Error: params.genome_list_composite must be a list containing the path of two fasta files"
    }

    ch_reads.view()
    ch_versions.view()

    /*
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
    */


}