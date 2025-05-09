process COMBINE_GENOMES {
    tag "$meta.id"
    label 'process_single'
    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(reference_fastas)

    output:
    tuple val(meta), path ("composite_genome.fa")   , emit: composite_genome
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat \\
        ${args} \\
        ${reference_fastas} > composite_genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat (GNU coreutils): \$(cat --version |& sed '1!d' | grep -E '[0-9].*')
    END_VERSIONS
    """

    stub:
    """
    touch composite_genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat (GNU coreutils): \$(cat --version |& sed '1!d' | grep -E '[0-9].*')
    END_VERSIONS
    """
}
