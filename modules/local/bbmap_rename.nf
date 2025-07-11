process BBMAP_RENAME {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01 bioconda::samtools=1.16.1 pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' :
        'biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' }"

    input:
    tuple val(meta), path(asm, stageAs: "input/*")
    val pref
    val minSize 

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: rename
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = pref ? "${pref}_${meta.id}_node_" : "${meta.id}_node_"

    """
    bbrename.sh in=$asm out=${meta.id}.fasta prefix=${prefix} minscaf=${minSize}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}