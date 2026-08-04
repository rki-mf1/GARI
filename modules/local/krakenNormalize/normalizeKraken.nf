process KRAKEN_NORMALIZE {
    label 'process_single'
    
    conda "conda-forge::python=3.9.21 conda-forge::ete3=3.1.3 conda-forge::pandas=2.2.3"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://caromanesco/gari/python_ete3_pandas:v.1.0.0' :
        'library://caromanesco/gari/python_ete3_pandas:v.1.0.0' }"

    input: 
        tuple val(metas), path(kraken)
        file(thresholds)
        path(ete3DB)

    output:
        tuple val(metas), path("*.classifiedreads.normalized.txt"), emit: report_norm
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def ete3_dbPath = ete3DB ? "-d $ete3DB" : ""
    def header = "sample,species,kraken2"
    def meta_list = metas instanceof List ? metas : [metas]
    def kraken_list = kraken instanceof List ? kraken : [kraken]
    def rows = meta_list.indices.collect { i -> "${meta_list[i].id},${meta_list[i].species},${kraken_list[i]}" }.join('\n')
    """
    echo "${header}\n${rows}" > samplesheet.csv

    krakenNorm.py  \\
        -t $thresholds \\
        -s samplesheet.csv \\
        $ete3_dbPath \\

    echo '"${task.process}":' > versions.yml
    echo "    python: \$(python --version | sed 's/Python //g')" >> versions.yml
    """
}