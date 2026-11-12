process QUARTO_REPORT {

    publishDir params.outdir, mode: 'copy'

    conda "environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/quarto_r-base_r-plotly_r-dt_pruned:c716339754baed25' :
        'community.wave.seqera.io/library/quarto_r-base_r-plotly_r-dt_pruned:c716339754baed25' }"

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
    # Set environment variables needed for Quarto rendering
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"

    # Without a _quarto.yml here, Quarto searches parent directories for a project and can
    # pick one outside the task (https://github.com/quarto-dev/quarto-cli/issues/14980)
    [ -e _quarto.yml ] || printf 'project:\\n  type: default\\n' > _quarto.yml

    # Fix Quarto for Apptainer (see https://community.seqera.io/t/confusion-over-why-a-tool-works-in-docker-but-fails-in-singularity-when-the-installation-doesnt-differ-i-e-using-wave-micromamba/1244)
    ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
    set +u
    if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
        source "\${ENV_QUARTO}"
    fi
    set -u

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