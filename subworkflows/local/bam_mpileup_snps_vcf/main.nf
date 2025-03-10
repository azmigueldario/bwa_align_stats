
include { TABIX_BGZIP           } from '../../../modules/nf-core/tabix/bgzip/main'
include { SAMTOOLS_FAIDX        } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_COVERAGE     } from '../../../modules/nf-core/samtools/coverage/main'  
include { SAMTOOLS_MPILEUP      } from '../../../modules/nf-core/samtools/mpileup/main'
include { BCFTOOLS_MPILEUP      } from '../../../modules/nf-core/bcftools/mpileup/main'      

workflow BAM_MPILEUP_SNPS_VCF {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_bai                  // channel: [ val(meta), [ bai ] ]
    ch_ref_fasta_gz         // channel: [ val(meta), [fasta] ]
    save_mpileup            // boolean: save mpileup file or not


    main:

    ch_versions = Channel.empty()

    //
    // Create inputs
    //

    TABIX_BGZIP(ch_ref_fasta_gz) 
    ch_ref_fasta = TABIX_BGZIP.out.output                               // channel with unzipped reference fasta  
    
    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output, [[], []] )                   // channel with indexed ref_fasta (no indexing intervals)
    
    ch_bam.join(ch_bai).set{ch_bam_bai}                                 // channel with bam and bai files together
    ch_mpileup_input = ch_bam.map{ meta, bam -> [ meta, bam, [] ] }     // channel for mpileup (no indexing intervals)

    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)             // capture versions
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // Calculate alignment coverage
    //
    
    if (!params.skip_coverage) {
        SAMTOOLS_COVERAGE(
            ch_bam_bai,
            ch_ref_fasta,
            SAMTOOLS_FAIDX.out.fai )
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)
    }
    
    //
    // Produce pileup document

        BCFTOOLS_MPILEUP( 
            ch_mpileup_input,
            ch_ref_fasta,
            save_mpileup)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    //
    // Call variants
    //

    emit:
    faidx       = SAMTOOLS_FAIDX.out.fai                // channel: [ val(meta), [ faidx ] ] 
    coverage    = SAMTOOLS_COVERAGE.out.coverage        // channel: [ val(meta), [ txt ] ] 
    pileup      = BCFTOOLS_MPILEUP.out.mpileup          // channel: [ val(meta), [ mpileup_vcf ] ]
    vcf         = BCFTOOLS_MPILEUP.out.vcf              // channel: [ val(meta), [ *vcf ] ] 
    vcf_index   = BCFTOOLS_MPILEUP.out.tbi              // channel: [ val(meta), [ tbi ] ]
    versions    = ch_versions                           // channel: [ versions.yml ]
}

