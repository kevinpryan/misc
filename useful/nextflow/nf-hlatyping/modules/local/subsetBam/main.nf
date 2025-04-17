process subsetBam{
    publishDir "$params.outdir/subsetBam"

    input:
    //tuple val(meta), path(bamfile), path(bamfileIndex)
    tuple val(meta), path(bam_bai) 
    path alt_chr6_contigs
    path hla_contigs
    path reference
    path fasta_bed
    path subset_regions
    output:
    tuple val(meta), path("*_subset.sorted.bam*"), emit: subsetbam
    path("*_subset.sorted.bam.flagstat")
    path("*_subset.sorted.bam.header") 
    shell:
    '''
    echo 'subsetting regions...'
    samtools view -h *.bam | \
    awk 'NR==FNR {a[$1]; next} /^@/ || $3 in a' !{subset_regions} - | \
    samtools view -Sb - > !{meta.sample}_subset1.bam 
    echo 'subsetting unmapped reads...'
    samtools view -o !{meta.sample}_subset2.bam -b !{meta.sample}_sorted_mdup.bam "*"
    echo 'merging bams...'
    samtools merge !{meta.sample}_subset.bam !{meta.sample}_subset1.bam !{meta.sample}_subset2.bam 
    echo 'sorting merged bams...'
    samtools sort -o !{meta.sample}_subset.sorted.bam !{meta.sample}_subset.bam
    echo 'indexing sorted merged bam...'
    sambamba index !{meta.sample}_subset.sorted.bam
    echo 'viewing header...'
    samtools view -H !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.header
    echo 'flagstat'
    samtools flagstat !{meta.sample}_subset.sorted.bam > !{meta.sample}_subset.sorted.bam.flagstat
    '''
}

//     samtools view -o ${meta.sample}_subset.bam -b *.bam chr6 chr6_ chrUn chrUn_ HLA "*"

// samtools view -o ${meta.sample}_subset.bam -b *.bam chr6:28509970-33480727 chr6* HLA* chrUn* "*"

/*
    samtools view -o ${meta.sample}_subset.bam -b ${bamfile} chr6 chr6_* HLA* "*"


    samtools view -o ${meta.sample}_subset.bam -b ${bamfile} -F 256 chr6:29500000-33200000 chr6_* HLA*

samtools view -o ${meta.sample}_subset1.bam -b ${bamfile} -L ${fasta_bed}
    samtools view -o ${meta.sample}_subset2.bam -b ${bamfile} "*"
    samtools view -o ${meta.sample}_subset3.bam -b ${bamfile} 'chr6:28509970-33480727'
    samtools merge ${meta.sample}_subset.bam ${meta.sample}_subset1.bam ${meta.sample}_subset2.bam ${meta.sample}_subset3.bam
    samtools sort -o ${meta.sample}_subset.sorted.bam ${meta.sample}_subset.bam
    sambamba index ${meta.sample}_subset.sorted.bam
    samtools view -H ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.header
    samtools flagstat ${meta.sample}_subset.sorted.bam > ${meta.sample}_subset.sorted.bam.flagstat
*/
