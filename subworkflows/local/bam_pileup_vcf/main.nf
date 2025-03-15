
include { BCFTOOLS_MPILEUP      } from '../../../modules/nf-core/bcftools/mpileup/main'      
include { FREEBAYES             } from '../../../modules/nf-core/freebayes/main'

workflow BAM_PILEUP_VCF {

    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai] ]
    ch_ref_fasta        // channel: [ val(meta), [ genome.fa ] ]
    ch_ref_fasta_fai    // channel: [ val(meta), [ genome.fasta.fai ] ]
    save_mpileup        // boolean: params.save_mpileup

    main:

    ch_versions = Channel.empty()

    // Variant calling and mpileup
    BCFTOOLS_MPILEUP( 
        ch_bam.map{meta, bam -> [ meta, bam, [] ] },
        ch_ref_fasta,
        save_mpileup)

    if (!params.skip_freebayes){
        FREEBAYES(
            ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, [], [], []] },
            ch_ref_fasta,
            ch_ref_fasta_fai,
            [[id:'null'], []],      // No specific list of samples to include
            [[id:'null'], []],      // No subpopulations specified
            [[id:'null'], []]       // Required for CNV analysis
        )
    }

    // Capture versions
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    emit:
    mpileup         = BCFTOOLS_MPILEUP.out.mpileup    // channel: [ val(meta), [*mpileup.gz] ]
    vcf_bcftools    = BCFTOOLS_MPILEUP.out.vcf        // channel: [ val(meta), [*vcf] ]
    vcf_freebayes   = FREEBAYES.out.vcf               // channel: [ val(meta), [*vcf] ]
    versions        = ch_versions                     // channel: [ versions.yml ]
}

