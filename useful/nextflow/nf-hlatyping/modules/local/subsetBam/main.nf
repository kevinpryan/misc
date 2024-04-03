process subsetBam{
    publishDir "$params.outdir/subsetBam", mode: 'copy'

    input:
    tuple val(meta), path(bamfile), path(bamfileIndex)
    path alt_chr6_contigs
    path hla_contigs
    path reference
    path fasta_bed
    output:
    tuple val(meta), path("*_subset.sorted.bam*"), emit: subsetbam
    path("*_subset.sorted.bam.flagstat")
    path("*_subset.sorted.bam.header") 
    script:
    """
    samtools view -o ${meta.sample}_subset1.bam -b ${bamfile} -L ${fasta_bed}
    samtools view -o ${meta.sample}_subset2.bam -b ${bamfile} "*"
    samtools view -o ${meta.sample}_subset3.bam -b ${bamfile} 'chr6:28509970-33480727'
    samtools merge ${meta.sample}_subset.bam ${meta.sample}_subset1.bam ${meta.sample}_subset2.bam ${meta.sample}_subset3.bam
    samtools sort -o ${meta.sample}_subset.sorted.bam ${meta.sample}_subset.bam
    sambamba index ${meta.sample}_subset.sorted.bam
    samtools view -H ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.header
    samtools flagstat ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.flagstat
    """
}
//     tuple val(meta), path("*_subset.sorted.bam"), path("*_subset.sorted.bam.bai"), emit: subsetbam
