process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01 bioconda::samtools=1.16.1 pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' :
        'biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' }"

    input:
    tuple val(meta), path(fastq), path(ref)

    output:
    tuple val(meta), path("*.covstats"), emit: covstats
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: log 
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"


    """
    bbmap.sh \\
        $input \\
        out=${prefix}.bam \\
        ref=$ref \\
        covstats=${prefix}.covstats \\
        $args \\
        -Xmx2g \\
        threads=$task.cpus

    cp .command.err ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
