
include { FASTP                                     } from '../../../modules/nf-core/fastp'
include { UNZIP as UNZIP_FASTA                      } from '../../../modules/nf-core/unzip/main'  
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx/main'
include { BWAMEM2_INDEX                             } from '../../../modules/nf-core/bwamem2/index/main'
include { GATK4_CREATESEQUENCEDICTIONARY            } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'

workflow PREPARE_INPUTS {

    take:
    ch_fastq_reads              // tuple [meta, [fastq1, fastq2] ] from '--input samplesheet.csv'
    ch_refgenome                // tuple [meta, fasta] from '--ref_genome path/to/fasta'
    ref_genome_path             // params.ref_genome
    fastp_save_trimmed_fail     // params.fastp_save_trimmed_fail
    fastp_save_merged           // params.fastp_save_merged

    main:

    ch_versions = Channel.empty()

    // Trim and process fastq files
    FASTP(
        ch_fastq_reads,
        [],         // No adapter fastas provided
        false,      // Keep reads passing trimming threshold
        fastp_save_trimmed_fail,
        fastp_save_merged    )

    // Unzip fasta if necessary
    if (ref_genome_path.endsWith(".gz")){
        UNZIP_FASTA(ch_refgenome)
        UNZIP_FASTA.out.unzipped_archive
            .map{ meta, dir ->
                [ meta, file("${dir}/*.{fna,fasta,fa}") ] }
            .set{ ch_ref_fasta }
        ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
    } else {
        ch_ref_fasta = ch_refgenome
    }

    // Create genome.fasta.fai and alignment index
    SAMTOOLS_FAIDX(ch_ref_fasta, [ [ id:'no_fai' ], [] ] )
    BWAMEM2_INDEX(ch_ref_fasta)

    // Create sequence dictionary for GATK4 processes
    GATK4_CREATESEQUENCEDICTIONARY(ch_ref_fasta)

    // Gather tool versions
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    reads           = FASTP.out.reads                           // channel: [ val(meta), [ fastq1, fastq2 ] ]
    ref_fasta       = ch_ref_fasta                              // channel: [ val(meta), [ genome.fa ] ]
    fasta_fai       = SAMTOOLS_FAIDX.out.fai                    // channel: [ val(meta), [ genome.fasta.fai ] ]
    bwamem2_index   = BWAMEM2_INDEX.out.index                   // channel: [ val(meta), [ index_dir ] ]
    dictionary      = GATK4_CREATESEQUENCEDICTIONARY.out.dict   // channel: [ val(meta), [*.dict] ]
    versions        = ch_versions                               // channel: [ versions.yml ]
}

