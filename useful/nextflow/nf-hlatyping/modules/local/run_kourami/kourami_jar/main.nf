
process RUN_KOURAMI_JAR{
    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(bam_bai)
    path kourami_graph
    output:
    tuple val(meta), path("*.result"), emit: kourami_result
    script:
    """
    java -jar /opt/wtsi-cgp/java/Kourami.jar.sh # ${meta.sample} "*.bam" -d ${kourami_panel} -r ${reference}
    """
}
//     HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}

