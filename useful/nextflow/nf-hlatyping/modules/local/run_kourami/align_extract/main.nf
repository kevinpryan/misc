
process RUN_KOURAMI_ALIGN_EXTRACT{
    publishDir "${params.outdir}/kourami/${meta.sample}", mode: 'copy'
    input:
    tuple val(meta), path(bam_bai)
    path kourami_panel
    path reference
    output:
    tuple val(meta), path("*on_KouramiPanel.bam), emit: kourami_alignment
    script:
    """
    bash alignAndExtract_hs38DH.sh ${meta.sample} "*.bam" -d ${kourami_panel} -r ${reference}
    """
}
//     HLA-LA.pl --BAM *.bam --graph /usr/local/bin/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT --workingDir . --sampleID ${meta.sample} --maxThreads ${task.cpus}

