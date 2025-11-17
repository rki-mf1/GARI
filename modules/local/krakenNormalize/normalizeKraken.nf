process KRAKEN_NORMALIZE {
    tag "$meta.id"
    label 'process_single'
    
    conda "conda-forge::python=3.9.21 conda-forge::ete3=3.1.3 conda-forge::pandas=2.2.3"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ete3%3A3.1.2' :
        'https://depot.galaxyproject.org/singularity/ete3%3A3.1.2' }"

    input: 
        tuple val(meta), path(kraken)
        file(thresholds)

    output:
        tuple val(meta), path("${meta.id}.classifiedreads.normalized.txt"), emit: report_norm

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    krakenNorm.py  \\
        -k $kraken \\
        -t $thresholds \\
        -s "${meta.species}" \\
        -o ${meta.id}.classifiedreads.normalized.txt
    """
    
}