process RUN_OPTITYPE{
    publishDir "optitype_out", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.sample}"), emit: optitype_call
    script:
    """
    OptiTypePipeline.py --input *.1.fq *.2.fq --verbose -d --outdir ${meta.sample}
    """
}