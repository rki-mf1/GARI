process DOWNLOAD_KRAKENDB {
    label 'process_single'

    conda "conda-forge::conda-forge::curl"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/curl:7.80.0' :
        'biocontainers/curl:7.80.0' }"

    input:

    output:
      path "./babykraken", emit: krakenDB
      path "versions.yml", emit: versions

    script:
    """
      curl -L https://github.com/MDU-PHL/babykraken/blob/master/dist/babykraken.tar.gz?raw=true | tar xz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n 1| sed 's/curl //g')
    END_VERSIONS

    """
    
}