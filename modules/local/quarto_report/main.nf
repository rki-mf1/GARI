process QUARTO_REPORT {

    tag "${qc_tsv.baseName}"
    publishDir params.outdir, mode: 'copy'

    conda "environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://bene1107/bioinfo/r-quarto-env:0.3' :
        'library://bene1107/bioinfo/r-quarto-env:0.3' }"

    input:
    path qc_tsv
    path report_qmd
    val report_title
    val report_author

    output:
    path "${qc_tsv.baseName}.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    quarto render ${report_qmd} \
        -P qc_file:${qc_tsv} \
        -P report_title:"${report_title}" \
        -P report_author:"${report_author}" \
        --output ${qc_tsv.baseName}.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}