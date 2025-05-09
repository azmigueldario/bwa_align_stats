
include { BCFTOOLS_MPILEUP              } from '../../../modules/nf-core/bcftools/mpileup/main'      
include { FREEBAYES                     } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_SORT                 } from '../../../modules/nf-core/bcftools/sort'
include { GATK4_HAPLOTYPECALLER         } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_VARIANTFILTRATION       } from '../../../modules/nf-core/gatk4/variantfiltration/main'


workflow BAM_PILEUP_VCF {

    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai] ]
    ch_ref_fasta        // channel: [ val(meta), [ genome.fa ] ]
    ch_ref_fasta_fai    // channel: [ val(meta), [ genome.fasta.fai ] ]
    ch_dictionary       // channel: [ val(meta), [*.dict] ]
    save_mpileup        // boolean: params.save_mpileup

    main:

    ch_versions         = Channel.empty()
    ch_vcf_gatk4        = Channel.empty()
    ch_vcf_freebayes    = Channel.empty()
    ch_vcf_bcftools     = Channel.empty()

    // Variant calling and mpileup
    BCFTOOLS_MPILEUP( 
        ch_bam.map{meta, bam -> [ meta, bam, [] ] },
        ch_ref_fasta,
        save_mpileup)

    BCFTOOLS_MPILEUP.out.vcf
        .map{ meta, vcf -> [ meta + [ variantcaller:'bcftools'], vcf ]}
        .set{ch_vcf_bcftools}
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    if (!params.skip_freebayes){
        FREEBAYES(
            ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, [], [], []] },
            ch_ref_fasta,
            ch_ref_fasta_fai,
            [[id:'null'], []],      // No specific list of samples to include
            [[id:'null'], []],      // No subpopulations specified
            [[id:'null'], []]       // Required for CNV analysis
        )
        BCFTOOLS_SORT(FREEBAYES.out.vcf)
        FREEBAYES.out.vcf
            .map{ meta, vcf -> [ meta + [ variantcaller:'freebayes'], vcf ] }
            .set{ch_vcf_freebayes}

        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
    }

    if (!params.skip_haplotypecaller) {
        GATK4_HAPLOTYPECALLER(
            ch_bam_bai.map{meta, bam, bai -> [ meta, bam, bai, [], [] ] },
            ch_ref_fasta,
            ch_ref_fasta_fai,
            ch_dictionary,
            [[id:'null'], []],      // No SNP reference database
            [[id:'null'], []]       // No SNP reference database index
        )
        GATK4_VARIANTFILTRATION(
            GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true),
            ch_ref_fasta,
            ch_ref_fasta_fai,
            ch_dictionary,
            [[id:'null'], []]           // No genome index is necessary for uncompressed fasta
        )
        
        GATK4_VARIANTFILTRATION.out.vcf
            .map{ meta, vcf -> [ meta + [ variantcaller:'gatk4'], vcf ] } 
            .set{ch_vcf_gatk4}

        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)
    }   

    emit:
    mpileup         = BCFTOOLS_MPILEUP.out.mpileup      // channel: [ val(meta), [*mpileup.gz] ]
    vcf_bcftools    = ch_vcf_bcftools                   // channel: [ val(meta), [*vcf] ]
    vcf_freebayes   = ch_vcf_freebayes                  // channel: [ val(meta), [*vcf] ]
    vcf_gatk4       = ch_vcf_gatk4                      // channel: [ val(meta), [*vcf] ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}

