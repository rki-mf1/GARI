process KRAKEN_NORMALIZE {
    tag "$meta.id"
    label 'process_single'
    
    conda "conda-forge::python=3.9.21 conda-forge::ete3=3.1.3 conda-forge::pandas=2.2.3"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://caromanesco/gari/python_ete3_pandas:v.1.0.0' :
        'library://caromanesco/gari/python_ete3_pandas:v.1.0.0' }"

    input: 
        tuple val(meta), path(kraken)
        file(thresholds)

    output:
        tuple val(meta), path("${meta.id}.classifiedreads.normalized.txt"), emit: report_norm
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    krakenNorm.py  \\
        -k $kraken \\
        -t $thresholds \\
        -s "${meta.species}" \\
        -o ${meta.id}.classifiedreads.normalized.txt
 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}