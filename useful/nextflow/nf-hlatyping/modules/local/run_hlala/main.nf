process RUN_HLALA{
    publishDir "${params.outdir}/hlala/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*.txt"), emit: hlala_call
    script:
    """
    HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --sampleID ${meta.sample} --maxThreads ${task.cpus}
    """
}
