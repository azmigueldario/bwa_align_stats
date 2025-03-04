
include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main' 
include { SAMTOOLS_MPILEUP      } from '../../../modules/nf-core/samtools/mpileup/main'

workflow BAM_MPILEUP_SNPS_VCF {

    take:
    ch_bam          // channel: [ val(meta), [ bam ] ]
    ch_ref_genome   // channel: [ val(meta), [fasta] ]


    main:

    ch_versions = Channel.empty()

        // empty intervals file, reproduce cardinality
    ch_bam
        .map{ meta, bam -> [ meta, bam, [] ] }
        .set{ch_mpileup_input}
    
    TABIX_BGZIP(ch_ref_genome)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)

    SAMTOOLS_MPILEUP( 
        ch_mpileup_input,
        TABIX_BGZIP.out.output.map{ meta, fasta -> fasta } )
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions)


    emit:
    pileup      = SAMTOOLS_MPILEUP.out.mpileup          // channel: [ val(meta), [ mpileup ] ]
    versions    = ch_versions                           // channel: [ versions.yml ]
}

