// Adapted from the vcf QC subworkflow in nf-core sarek

include { BCFTOOLS_STATS                  } from '../../../modules/nf-core/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_QUAL       } from '../../../modules/nf-core/vcftools/main'

workflow VCF_STATS_REPORTS {

    take:
    vcf_gatk4           // All: channel [ val(meta), (*.vcf)]
    vcf_freebayes
    vcf_bcftools

    main:

    ch_versions = Channel.empty()

    ch_vcf_all = vcf_gatk4
        .mix(vcf_freebayes)
        .mix(vcf_bcftools)
    
    BCFTOOLS_STATS(
        ch_vcf_all.map{ meta, vcf -> [ meta, vcf, [] ] }, 
        [[:],[]],       // No target regions, targets, exons, samples, or ref_fasta
        [[:],[]], 
        [[:],[]], 
        [[:],[]], 
        [[:],[]])
    
    VCFTOOLS_QUAL(ch_vcf_all, [], [])
    VCFTOOLS_SUMMARY(ch_vcf_all, [], [])

    // capture versions
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_QUAL.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats              // channel: [ val(meta), (*stats.txt) ]
    vcftools_qual           = VCFTOOLS_QUAL.out.tstv_qual           // channel: [ val(meta), (*TsTv.qual) ]
    vcftools_summary        = VCFTOOLS_SUMMARY.out.filter_summary   // channel: [ val(meta), (*FILTER.sumamry) ]
    versions                = ch_versions                           // channel: [ versions.yml ]
}

