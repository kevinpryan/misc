process NORMALISE_VCFS {
    label 'bcftools_container'
    publishDir "$params.outdir/norm_vcf", mode: 'copy'
    input:
    tuple val(patient), path(vcf)
    path(fasta)
    output:
    tuple val(patient), path("*_normalised.vcf.gz"), emit: vcf_norm
    script:
    """
    bcftools norm --fasta-ref $fasta --output-type z --output ${vcf.simpleName}_normalised.vcf.gz $vcf
    """
}


process VCF2BED {
    label 'bedops_container'
    publishDir "$params.outdir/vcf2bed", mode: 'copy'
    input:
    tuple val(patient), path(vcf_norm)
    output:
    tuple val(patient), path("*.bed"), emit: bed_from_vcf
    script:
    """
    vcf2bed < <( gunzip -c $vcf_norm ) | cut -f 1-3 > ${patient}.bed
    """
}

process SUBSET_CRAM {
    label 'samtools_container'
    publishDir "$params.outdir/subset_cram/${patient}", mode: 'copy'
    input:
    tuple val(patient), path(cram), path(bed)
    path(fasta)
    output:
    tuple val(patient), path("*_subset_mutect2.bam"), path("*_subset_mutect2.bam.bai")
    script:
    """
    samtools view --reference $fasta -bL $bed -o ${patient}_subset_mutect2.bam $cram
    samtools index ${patient}_subset_mutect2.bam
    """
}
