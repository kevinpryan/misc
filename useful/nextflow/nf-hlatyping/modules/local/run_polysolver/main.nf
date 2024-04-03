process RUN_POLYSOLVER{
    publishDir "${params.outdir}/polysolver/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("winners.hla.nofreq.txt"), emit: polysolver_call
    tuple val(meta), path("counts*")
    path("check.status.out.txt")
    script:
    """
    /home/polysolver/scripts/shell_call_hla_type *.bam Unknown 0 hg38 ILMFQ 0 ./
    """
}

