process EXTRACT_MAPPING_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path (reference_fastas)

    output:
    tuple val(meta), path("*.tsv")      , emit: mapping_stats
    path("versions.yml")                  , emit: versions          

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_a = reference_fastas[0]
    def genome_b = reference_fastas[1]
    """
    # Get the length of each genome
    genome_a_length=\$(grep -v '^>' ${genome_a} | tr -d '\\n' | wc -c)
    genome_b_length=\$(grep -v '^>' ${genome_b} | tr -d '\\n' | wc -c)

    # Calculate total mapped reads
    total_mapped=\$(samtools view -F 0x904 -c ${bam})

    # Calculate reads mapped to each genome
    genome_a_mapped=\$(samtools view -F 0x904 ${args} ${bam} | awk -v a=\$genome_a_length '\$4 <= a' | wc -l)
    genome_b_mapped=\$(samtools view -F 0x904 ${args} ${bam} | awk -v a=\$genome_a_length '\$4 > a' | wc -l)

    # Calculate percentages
    genome_a_percent=\$(awk "BEGIN {printf \\"%.2f\\", (\$genome_a_mapped / \$total_mapped) * 100}")
    genome_b_percent=\$(awk "BEGIN {printf \\"%.2f\\", (\$genome_b_mapped / \$total_mapped) * 100}")

    # Write results to tab separated file
    echo -e "${prefix}\t\${total_mapped}\t\${genome_a_mapped}\t\${genome_b_mapped}\t\${genome_a_percent}\t\${genome_b_percent}" > ${prefix}_mapping_stats.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mapping_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
