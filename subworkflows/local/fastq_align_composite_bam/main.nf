include { COMBINE_GENOMES                   } from '../../../modules/local/combinegenomes/main.nf'
include { BWAMEM2_INDEX as COMPOSITE_INDEX  } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as COMPOSITE_ALIGN    } from '../../../modules/nf-core/bwamem2/mem/main'
include { EXTRACT_MAPPING_STATS             } from '../../../modules/local/extractmappingstats/main.nf'
include { CONCATENATE_MAPPING_STATS         } from '../../../modules/local/concatenatemappingstats/main.nf'


workflow FASTQ_ALIGN_COMPOSITE_BAM {

    take:
    ch_reads        // channel: [mandatory] meta, reads
    genome_list     // params.genome_list_composite
    sort_bam        // params.sort_bam

    main:

    ch_versions = Channel.empty()

    // Check input reference genomes and create input channel
    if (genome_list && genome_list instanceof List && genome_list.size() == 2 ) {

        if (!file(genome_list[0]).exists() || !file(genome_list[1]).exists() ) {
            error "Error: At least one of the paths specified in ${genome_list} does not exist"
        }
        Channel
            .fromList(genome_list)
            .collect()
            .map { [ [id:'input_genomes_composite'], [it[0], it[1]] ] }
            .set{input_composite_ch}
    }
    else { 
        error "Error: params.genome_list_composite must be a list containing valid paths to 2 fasta files"
    }

    // Combine and index reference genomes
    COMBINE_GENOMES(input_composite_ch)
    COMPOSITE_INDEX(COMBINE_GENOMES.out.composite_genome)
    COMPOSITE_ALIGN(
        ch_reads,
        COMPOSITE_INDEX.out.index,
        COMBINE_GENOMES.out.composite_genome,
        sort_bam
    )   

    // Summarize alignment stats
    EXTRACT_MAPPING_STATS(
        COMPOSITE_ALIGN.out.bam,
        input_composite_ch
    )
    EXTRACT_MAPPING_STATS.out.mapping_stats
        .map{ _meta, tsv -> [ [id:'mapping_stats'], tsv ] }
        .groupTuple()
        .collect()
        .set{mapping_reports_ch}
    CONCATENATE_MAPPING_STATS(
        mapping_reports_ch,
        input_composite_ch
    )
    
    // Versions
    ch_versions = ch_versions.mix(COMBINE_GENOMES.out.versions)
    ch_versions = ch_versions.mix(COMPOSITE_INDEX.out.versions)
    ch_versions = ch_versions.mix(COMPOSITE_ALIGN.out.versions)
    ch_versions = ch_versions.mix(EXTRACT_MAPPING_STATS.out.versions)
    ch_versions = ch_versions.mix(CONCATENATE_MAPPING_STATS.out.versions)

    emit:
    report   = CONCATENATE_MAPPING_STATS.out.summary        // channel: [val(meta), [summary.tsv]]
    versions = ch_versions                                  // channel: [ versions.yml ]
}