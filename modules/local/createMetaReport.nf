process STAT_SUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
      tuple val(meta), file(ASM), file(BUSCO), file(FASTANI), file(STATS), file(FASTP), file(KRAKEN), file(KRAKEN_ASM), file(BBMAP), file(CHECKM)
      file(refLens)
      file(thresholds)
      path(krakenDB)

    output:
      tuple val(meta), path("./${meta.id}_QC.json"), emit: json
      path "versions.yml", emit: versions

    script:
    """
    createReport.py \\
      --b $BUSCO \\
      --cm $CHECKM \\
      --f $FASTANI \\
      --s $STATS \\
      --p $FASTP \\
      --kr $KRAKEN \\
      --ka $KRAKEN_ASM \\
      --r ${refLens} \\
      --c ${BBMAP} \\
      --o ${meta.id}_QC.json \\
      --sp "${meta.species}" \\
      --gv ${workflow.manifest.version} \\
      --ga ${params.assembler} \\
      --gfDB ${params.reference} \\
      --gbDB ${params.busco_lin} \\
      --gkDB ${krakenDB} \\
      --th ${thresholds} \\
      --ck ${params.classify_kraken}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
