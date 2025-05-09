
process CONCATENATE_MAPPING_STATS {
    tag "$meta.id"
    label 'process_single'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(reports)
    tuple val(meta2), path(reference_fastas)


    output:
    tuple val(meta), path("*summary.tsv")      , emit: summary
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_a = reference_fastas[0].getSimpleName()
    def genome_b = reference_fastas[1].getSimpleName()
    """
    # concatenate all results
    echo -e ${args} "sample_id\ttotal_reads_mapped\t${genome_a}_mapped\t${genome_b}_mapped\t${genome_a}_percent\t${genome_b}_percent" > temp.tsv
    cat *_mapping_stats.tsv >> temp.tsv

    # new column with name of genome where the majority of reads mapped
    awk -v var_a="${genome_a}" \
        -v var_b="${genome_b}" '
        {
        if (NR == 1) 
            print \$0,"likely_taxa";
        else  if (\$3>=\$4) 
            print \$0,var_a;
        else 
            print \$0,var_b 
        }
    ' temp.tsv > alignment_summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat (GNU coreutils): \$(cat --version |& sed '1!d' | grep -E '[0-9].*')
    END_VERSIONS
    """

    stub:
    """
    touch alignment_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat (GNU coreutils): \$(cat --version |& sed '1!d' | grep -E '[0-9].*')
    END_VERSIONS
    """
}
