process RUN_OPTITYPE{
    publishDir "${params.outdir}/optitype/${meta.sample}", mode: 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/optitype:1.3.5--0' :
        'quay.io/biocontainers/optitype:1.3.5--0' }"
    input:
    tuple val(meta), path(reads)
    //val dna_rna
    output:
    tuple val(meta), path("*.tsv"), emit: optitype_call
    tuple val(meta), path("*.pdf")
    script:
    """
    OptiTypePipeline.py --input *.1.fq *.2.fq --verbose --${meta.seq_type} --outdir ${meta.sample}
    cp "${meta.sample}"/*/*.tsv .
    cp "${meta.sample}"/*/*.pdf .
    """
}
