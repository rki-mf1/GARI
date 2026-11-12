process CREATE_REPORT {
    label 'process_single'

    conda "conda-forge::pandas=1.2.5 conda-forge::openpyxl=3.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pandas_openpyxl:acd440c10fea94e1' :
        'community.wave.seqera.io/library/pandas_openpyxl:acd440c10fea94e1' }"

    input:
      file ('*')
      val prefix

    output:
      path '*.tsv'        , emit: tsv
      path '*.xlsx'       , emit: xlxs
      path "versions.yml" , emit: versions

    script:
    """
    summarizeReports.py \\
        --g . \\
        --p $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
