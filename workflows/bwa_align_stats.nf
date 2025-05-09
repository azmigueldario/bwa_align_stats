/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_bwa_align_stats_pipeline'
include { PREPARE_INPUTS            } from '../subworkflows/local/prepare_inputs'
include { FASTQ_ALIGN_COMPOSITE_BAM } from '../subworkflows/local/fastq_align_composite_bam/main.nf'
include { FASTQ_ALIGN_DNA           } from '../subworkflows/nf-core/fastq_align_dna'
include { BAM_PROCESSING_QC_STATS   } from '../subworkflows/local/bam_processing_qc_stats'  
include { BAM_PILEUP_VCF            } from '../subworkflows/local/bam_pileup_vcf/main.nf'
include { VCF_STATS_REPORTS         } from '../subworkflows/local/vcf_stats_reports/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BWA_ALIGN_STATS {

    take:
    ch_fastq_reads              // tuple [meta, [fastq1, fastq2] ] from '--input samplesheet.csv'
    ch_refgenome                // tuple [meta, fasta] from '--ref_genome fasta'
    ref_genome_path             // params.ref_genome
    fastp_save_trimmed_fail     // params.fastp_save_trimmed_fail
    fastp_save_merged           // params.fastp_save_merged
    composite_references_input  // params.genome_list_composite
    aligner                     // params.aligner
    sort_bam                    // params.sort_bam
    save_mpileup                // params.save_mpileup

    main:

    ch_versions = Channel.empty()

    // Add readgroups to the metamap
    input_fastq = ch_fastq_reads.map{meta, files -> addReadgroupToMeta(meta, files)}

    // 
    // FastQ processing, prepare inputs
    //

    PREPARE_INPUTS(
        input_fastq,
        ch_refgenome,
        ref_genome_path,
        fastp_save_trimmed_fail,
        fastp_save_merged
    )
    ch_versions = ch_versions.mix(PREPARE_INPUTS.out.versions)

    //
    // Classify assemblage of input reads
    //
    FASTQ_ALIGN_COMPOSITE_BAM(
        PREPARE_INPUTS.out.reads,
        composite_references_input,
        sort_bam
    )
    
    //
    // Align reads to reference genome
    //

    FASTQ_ALIGN_DNA(
        PREPARE_INPUTS.out.reads,
        PREPARE_INPUTS.out.bwamem2_index,
        PREPARE_INPUTS.out.ref_fasta,
        aligner,
        sort_bam )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)

    // 
    // Processing alignment file
    //

    BAM_PROCESSING_QC_STATS(
        FASTQ_ALIGN_DNA.out.bam,
        PREPARE_INPUTS.out.ref_fasta,
        PREPARE_INPUTS.out.fasta_fai,
        PREPARE_INPUTS.out.dictionary
    )
    ch_versions = ch_versions.mix(BAM_PROCESSING_QC_STATS.out.versions)

    //
    // Variant calling and mpileup
    //

    BAM_PILEUP_VCF(
        BAM_PROCESSING_QC_STATS.out.bam,
        BAM_PROCESSING_QC_STATS.out.bam_bai,
        PREPARE_INPUTS.out.ref_fasta,
        PREPARE_INPUTS.out.fasta_fai,
        PREPARE_INPUTS.out.dictionary,
        save_mpileup
    )
    ch_versions = ch_versions.mix(BAM_PILEUP_VCF.out.versions)

    VCF_STATS_REPORTS(
        BAM_PILEUP_VCF.out.vcf_gatk4,
        BAM_PILEUP_VCF.out.vcf_bcftools,
        BAM_PILEUP_VCF.out.vcf_freebayes
    )
    ch_versions = ch_versions.mix(VCF_STATS_REPORTS.out.versions)

    //
    // Collate and save software versions
    //

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'bwa_align_stats_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    

    emit:
    index          = PREPARE_INPUTS.out.bwamem2_index       // channel: [ val(meta), path(bwamem2_dir)]
    // reports        = BAM_PROCESSING_QC_STATS.out.reports    // channel: [ val(meta), [report1, report2 ...]]
    mpileup        = BAM_PILEUP_VCF.out.mpileup            // channel: [ val(meta), path(vcf_mpileup)]
    vcf_freebayes  = BAM_PILEUP_VCF.out.vcf_freebayes      // channel: [ val(meta), path(vcf_mpileup)]
    versions       = ch_collated_versions                  // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTION FOR READ GROUPS (modified from nf-core/sarek)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def flowcell = flowcellLaneFromFastq(files[0])

    // Check if flowcell ID matches
    if ( flowcell && flowcell != flowcellLaneFromFastq(files[1]) ){
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def sample_lane_id = flowcell ? "${flowcell}.${meta.id}" : "${meta.id}"

    // Don't use a random element for ID, it breaks resuming
    def read_group = "\"@RG\\tID:${sample_lane_id}\\t${CN}PU:${sample_lane_id}\\tSM:${meta.id}\\tLB:${meta.id}\\tDS:${meta.ref_genome_version}\\tPL:ILLUMINA\""
    meta  = meta - meta.subMap('lane') + [read_group: read_group.toString()]
    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // First line of FASTQ file contains sequence identifier plus optional description
    def firstLine = readFirstLineOfFastq(path)
    def flowcell_id = null

    // Expected format from ILLUMINA
    // cf https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    // Five fields:
    // @<instrument>:<lane>:<tile>:<x-pos>:<y-pos>...
    // Seven fields or more (from CASAVA 1.8+):
    // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>..."

    def fields = firstLine ? firstLine.split(':') : []
        if (fields.size() == 5) {
            // Get the instrument name as flowcell ID
            flowcell_id = fields[0].substring(1)
        } else if (fields.size() >= 7) {
            // Get the actual flowcell ID
            flowcell_id = fields[2]
        } else if (fields.size() != 0) {
            log.warn "FASTQ file(${path}): Cannot extract flowcell ID from ${firstLine}"
        }
        return flowcell_id
}

// Get first line of a FASTQ file
def readFirstLineOfFastq(path) {
    def line = null
    try {
        path.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            line = buffered.readLine()
            assert line.startsWith('@')
        }
    } catch (Exception e) {
        log.warn "FASTQ file(${path}): Error streaming"
        log.warn "${e.message}"
    }
    return line
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
