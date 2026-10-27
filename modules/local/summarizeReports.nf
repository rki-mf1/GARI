process CREATE_REPORT {
    label 'process_single'

    conda "conda-forge::pandas=1.2.5 conda-forge::openpyxl"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://mdriller/gari/pandas_openpyxl:latest' :
        'library://mdriller/gari/pandas_openpyxl:latest' }"

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
