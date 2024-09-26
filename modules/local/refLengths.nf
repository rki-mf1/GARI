process REF_LENGTHS {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
        path(refList)

    output:
        path("./refLengths.tsv"), emit: table
        path("./refList_firstCol.tsv"), emit: refListfastANI
        path "versions.yml", emit: versions

    script:
    """
    cut -f1 $refList > refList_firstCol.tsv

    createRefLenTable.py --r refList_firstCol.tsv --o refLengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
