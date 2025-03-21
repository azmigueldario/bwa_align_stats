
include { BCFTOOLS_MPILEUP              } from '../../../modules/nf-core/bcftools/mpileup/main'      
include { FREEBAYES                     } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_SORT                 } from '../../../modules/nf-core/bcftools/sort'
include { GATK4_HAPLOTYPECALLER         } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow BAM_PILEUP_VCF {

    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    ch_bam_bai          // channel: [ val(meta), [ bam ], [bai] ]
    ch_ref_fasta        // channel: [ val(meta), [ genome.fa ] ]
    ch_ref_fasta_fai    // channel: [ val(meta), [ genome.fasta.fai ] ]
    ch_dictionary       // channel: [ val(meta), [*.dict] ]
    save_mpileup        // boolean: params.save_mpileup

    main:

    ch_versions = Channel.empty()

    // Variant calling and mpileup
    BCFTOOLS_MPILEUP( 
        ch_bam.map{meta, bam -> [ meta, bam, [] ] },
        ch_ref_fasta,
        save_mpileup)
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
            /*
                cnn_ch = GATK4_HAPLOTYPECALLER.out.vcf
                    .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
                    .map{ meta, vcf, tbi -> [ meta, vcf, tbi, [], [] ]}
                GATK4_CNNSCOREVARIANTS(
                    cnn_ch,
                    ch_ref_fasta.map{ meta, fasta -> fasta},
                    ch_ref_fasta_fai.map{ meta, fai -> fai},
                    ch_dictionary.map{meta, dict -> dict},
                    [],
                    []
                )
                filter_ch = GATK4_CNNSCOREVARIANTS.out.vcf
                    .join(GATK4_CNNSCOREVARIANTS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
                    .map{ meta, vcf, tbi -> [ meta, vcf, tbi, [] ] }
                GATK4_FILTERVARIANTTRANCHES(
                    filter_ch,
                    [],                     // No SNP reference database
                    [],                     // No SNP reference database index
                    ch_ref_fasta.map{ meta, fasta -> fasta},
                    ch_ref_fasta_fai.map{ meta, fai -> fai},
                    ch_dictionary.map{meta, dict -> dict}
                )
                if(params.joint_variant_calling){
                    
                    //  Create input channel for merging 
                    gvcf_tbi_intervals = GATK4_HAPLOTYPECALLER.out.vcf
                        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
                        .map{ 
                            meta, vcf, tbi ->
                            [ [id: 'joint_variant_calling', intervals:'no_intervals'], vcf, tbi ] }       
                        .groupTuple(by:0)                                                               // Merge all .vcf together
                        .map{ meta, vcf, tbi -> [meta, vcf, tbi, [], [], []] }                           // No interval_file, 0 as num_intervals

                    GATK4_GENOMICSDBIMPORT(gvcf_tbi_intervals, false, false, false)
                    GATK4_GENOTYPEGVCFS(
                        GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{ meta, db -> [ meta, db, [], [], [] ]},
                        ch_ref_fasta,
                        ch_ref_fasta_fai,
                        ch_dictionary,
                        [[id:'null'], []],      // No SNP reference database
                        [[id:'null'], []]       // No SNP reference database index
                    )
                }
                */
    }

    emit:
    mpileup         = BCFTOOLS_MPILEUP.out.mpileup    // channel: [ val(meta), [*mpileup.gz] ]
    vcf_bcftools    = BCFTOOLS_MPILEUP.out.vcf        // channel: [ val(meta), [*vcf] ]
    vcf_freebayes   = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), [*vcf] ]
    vcf_gatk4       = GATK4_HAPLOTYPECALLER.out.vcf   // channel: [ val(meta), [*vcf] ]
    versions        = ch_versions                     // channel: [ versions.yml ]
}

