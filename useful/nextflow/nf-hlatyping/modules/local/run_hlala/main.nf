process RUN_HLALA{
    publishDir "${params.outdir}/hlala/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path graphdir
    output:
    tuple val(meta), path("${meta.sample}/hla/R1_bestguess_G.txt"), emit: hlala_call
    script:
    """
    HLA-LA.pl --BAM *.bam --customGraphDir ${graphdir} --graph PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}
    """
}
//     HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}

