process RUN_OPTITYPE{
    publishDir "${params.outdir}/optitype/${meta.sample}", mode: 'copy'
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
