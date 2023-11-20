#!/usr/bin/env nextflow
params.outdir = "."

process FILTER_VCF {
    publishDir "$params.outdir/vcf_filter_snps_only", mode: 'copy'
    input:
    tuple val(id), path(vcf)
    output:
    path("*filtered.vcf"), emit: ch_vcf
    name = vcf.getSimpleName()
    script:
    """
    vcftools --vcf $vcf --remove-indels --recode --recode-INFO-all --out ${name}.intersect.snps_only.vcf
    """
}

workflow {
    Channel.fromPath("samplesheet.csv")
    | splitCsv( header: true )
    | map { row ->
        [row.id, [file(row.vcf)]]
    }
    | FILTER_VCF
}
