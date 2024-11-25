process SKESA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::skesa=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa%3A2.5.1--hdcf5f25_1' :
        'bbiocontainers/skesa:2.5.1--hdcf5f25_1' }"

    input:
    tuple val(meta), path(illumina)

    output:
    tuple val(meta), path("*.fa")      , emit: contigs
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def reads = illumina ? ( meta.single_end ? "--reads $illumina" : "--reads ${illumina[0]},${illumina[1]}" ) : ""

    """
    skesa \\
        $args \\
        --cores $task.cpus \\
        --memory $maxmem \\
        $reads \\
        --contigs_out ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skesa: 2.5.1
    END_VERSIONS
    """
}
