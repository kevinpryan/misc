#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//workflow {
    
//ch_mutect = Channel.fromFilePairs("/home/kevin/Documents/PhD/wes_bc/nfcore_sarek_results/annotation/mutect2/**.ann.vcf.gz*",  checkIfExists: true)
//ch_mutect.view {"value: $it"} 
//}

//ch_mutect = Channel.fromPath("/home/kevin/Documents/PhD/wes_bc/nfcore_sarek_results/annotation/mutect2/**.ann.vcf.gz",  checkIfExists: true)
//ch_mutect.view()

process BCFTOOLS_FILTER_MUTECT2 {
    publishDir "$params.outdir/mutect2_filter", mode: 'copy'
    input:
    path(mutect) // ch_hlatyping

    output:
    path("*.vcf.gz"), emit: ch_mutect_filter

    script:
    prefix = mutect.simpleName 
    """
    bcftools view -f PASS ${mutect} > ${prefix}.mutect2.filtered_snpEff_VEP.ann.only.PASS.vcf.gz
    """
}

process MUTECT2_EXTRACT_MISSENSE {
    publishDir "$params.outdir/mutect2_filter_missense", mode: 'copy'
    input:
    path(mutect) // ch_hlatyping

    output:
    path("*.vcf.gz"), emit: ch_mutect_keep_missense

    script:
    prefix = mutect.simpleName 
    """
    zgrep "missense_variant" ${mutect} > ${prefix}.filtered_snpEff_VEP.ann.pass.missense.vcf.gz
    """
}
// ch_mutect_pair = Channel.fromPath("/home/kevin/Documents/PhD/wes_bc/nfcore_sarek_results/annotation/mutect2/**{.ann.vcf.gz,.ann.vcf.gz.tbi}",  checkIfExists: true)

ch_mutect = Channel.fromPath("/home/kevin/Documents/PhD/wes_bc/nfcore_sarek_results/annotation/mutect2/**.ann.vcf.gz",  checkIfExists: true)

workflow {
    BCFTOOLS_FILTER_MUTECT2(ch_mutect)
    MUTECT2_EXTRACT_MISSENSE(BCFTOOLS_FILTER_MUTECT2.out.ch_mutect_filter)
}

