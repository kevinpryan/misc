#!/usr/bin/env nextflow

process fasta_index_bed{
    publishDir "$params.outdir/fasta_bed"

    input:
    path reference
    val chromosome
    val reference_basename
    output:
    path("*.bed"), emit: fasta_bed

    script:
    """
    samtools faidx ${reference_basename}.fa
    awk 'BEGIN {FS="\\t"}; {print \$1 FS "0" FS \$2}' ${reference_basename}.fa.fai > ${reference_basename}.full.fa.bed
    grep ${chromosome} ${reference_basename}.full.fa.bed > ${reference_basename}.fa.bed
    """
}
process bwa_mem_align_alt{
    publishDir "$params.outdir/align"

    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename

    output:
    tuple val(meta), path("*.bwamem.sam"), emit: samfile
    path("*.bwamem.sam.header")
    path("*.bwamem.sam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    samtools view -H ${meta.sample}.bwamem.sam > ${meta.sample}.bwamem.sam.header
    samtools flagstat ${meta.sample}.bwamem.sam > ${meta.sample}.bwamem.sam.flagstat
    """
}

process bwa_mem_postalt{
    publishDir "$params.outdir/postalt"

    input:
    path reference
    tuple val(meta), path(samfile)
    val reference_basename
    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.sam.header")
    path("*_postalt.sam.flagstat")
    path("*_postalt.bam.flagstat")

    script:
    """
    k8-linux /usr/local/bin/bwa-0.7.15/bwa-postalt.js \
    ${reference_basename}.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam
    samtools view -H ${meta.sample}_postalt.sam > ${meta.sample}_postalt.sam.header
    samtools flagstat ${meta.sample}_postalt.sam > ${meta.sample}_postalt.sam.flagstat
    samtools view -bh -o ${meta.sample}_postalt.bam ${meta.sample}_postalt.sam
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    """
}

process bwa_mem_align_alt_postalt{
    publishDir "$params.outdir/bwa-aln-postalt"
    input:
    path reference
    tuple val(meta), path(reads)
    val reference_basename
    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.bam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference_basename}.fa ${reads} > "${meta.sample}.bwamem.sam"
    k8-linux /usr/local/bin/bwa-0.7.15/bwa-postalt.js \
    ${reference_basename}.fa.alt \
    ${meta.sample}.bwamem.sam > ${meta.sample}_postalt.sam |\
    samtools view -bh -o ${meta.sample}_postalt.bam -
    rm "${meta.sample}.bwamem.sam"  "${meta.sample}_postalt.sam"
    samtools flagstat ${meta.sample}_postalt.bam > ${meta.sample}_postalt.bam.flagstat
    """
}
// ${meta.sample}_postalt.sam

process samtools_sort{
    publishDir "$params.outdir/sort"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sortedbam

    script:
    """
    samtools sort -o ${meta.sample}.sorted.bam ${bamfile}
    """
}

process samtools_index{
    publishDir "$params.outdir/index"

    input:
    tuple val(meta), path(sortedbam)

    output:
    tuple val(meta), path("*.sorted.bam.bai"), emit: bam_indexed

    script:
    """
    samtools index ${sortedbam}
    """
}   

process markduplicates{
    publishDir "$params.outdir/markduplicates"

    input:
    tuple val(meta), path(sortedbam), path(sortedbam_index)

    output:
    tuple val(meta), path("*_sorted_mdup.bam"), path("*_sorted_mdup.bam.bai"), emit: markdupbam
    path("*.sorted.mdup.bam.flagstat")

    script:
    """
    bammarkduplicates I=${sortedbam} O=${meta.sample}_sorted_mdup.bam index=1 rmdup=0
    samtools flagstat ${meta.sample}_sorted_mdup.bam > ${meta.sample}.sorted.mdup.bam.flagstat
    """
}

process extractContigs {
    publishDir "$params.outdir/extractContigs"

    input:
    path hlatypes
    path reference
    val reference_basename
    val chromosome
    output:
    path("bwakit-alt_contigs.txt"), emit: alt_contigs
    path("bwakit-hla_contigs.txt"), emit: hla_contigs
    script:
    """
    grep ${chromosome} ${reference_basename}.fa.alt | awk -F '\\t' '{print \$1}' > bwakit-alt_contigs.txt
    awk -F '\\t' '{print \$3}' ${hlatypes} > bwakit-hla_contigs.txt
    """
}

process extractContigsTest {
    publishDir "$params.outdir/extractContigs"

    input:
    path hlatypes
    path reference
    val reference_basename
    output:
    path("bwakit-alt_chr19_contigs.txt"), emit: alt_chr19_contigs
    path("bwakit-hla_contigs.txt"), emit: hla_contigs
    script:
    """
    grep 'chr19_' ${reference_basename}.fa.alt | awk -F '\\t' '{print \$1}' > bwakit-alt_chr19_contigs.txt
    awk -F '\\t' '{print \$3}' ${hlatypes} > bwakit-hla_contigs.txt
    """
}


process subsetBam{
    publishDir "$params.outdir/subsetBam"

    input:
    tuple val(meta), path(bamfile), path(bamfileIndex)
    path alt_chr6_contigs
    path hla_contigs
    path reference
    path fasta_bed
    output:
    tuple val(meta), path("*_subset.sorted.bam"), path("*_subset.sorted.bam.bai"), emit: subsetbam
    path("*_subset.sorted.bam.flagstat")
    path("*_subset.sorted.bam.flagstat") 
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

process subsetBamtest{
    publishDir "$params.outdir/subsetBam"

    input:
    tuple val(meta), path(bamfile), path(bamfileIndex)
    path alt_chr19_contigs
    path hla_contigs
    path reference
    path fasta_bed
    output:
    tuple val(meta), path("*_subset.sorted.bam"), path("*_subset.sorted.bam.bai"), emit: subsetbam
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
//     samtools view -o ${meta.sample}_subset.bam -M -b ${bamfile} '*' -L chr19_chr19_KI270866v1_alt.fa.bed

//    awk 'BEGIN {FS="\\t"}; {print \$1 FS "0" FS \$2}' chr19_chr19_KI270866v1_alt.fa.fai > chr19_chr19_KI270866v1_alt.fa.bed

//     samtools view -o ${meta.sample}_subset.bam -b ${bamfile} 'chr6:28509970-33480727' ${alt_chr19_contigs} ${hla_contigs} '*'

process bam2fastq{
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(subsetbam), path(subsetbam_bai)
    output:
    tuple val(meta), path("*.fq"), emit: subsetfastq

    script:
    """
    sambamba view \
    -f "bam" -h -p -l 0 -t ${task.cpus} \
    ${subsetbam} | sambamba sort -p -n -t ${task.cpus} -o - /dev/stdin |
    samtools fastq /dev/stdin \
    -1 "${meta.sample}_subset.1.fq" \
    -2 "${meta.sample}_subset.2.fq" \
    -0 /dev/null -s /dev/null -n
    """
}

process realignwithoutAlt{
    publishDir "$params.outdir/realignwithoutAlt"
    
    input:
    tuple val(meta), path(reads)
    path reference
    val reference_basename

    output:
    tuple val(meta), path("*_realign.sam"), emit: realignbam
    path("*_realign.sam.header")
    path("*_realign.sam.flagstat")

    script:
    """
    bwa mem -t ${task.cpus} -j ${reference_basename}.fa ${reads} > ${meta.sample}_realign.sam
    samtools view -H ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.header
    samtools flagstat ${meta.sample}_realign.sam > ${meta.sample}_realign.sam.flagstat
    """
}

process reheaderChr{
    publishDir "$params.outdir/reheaderChr"

    input:
    tuple val(meta), path(bamfile), path(bamfile_idx)

    output:
    //tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), path("*_reheader.bam.flagstat"), path("*_reheader.bam.header"), emit: reheaderbam
    tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), emit: reheaderbam
    path("*_reheader.bam.flagstat")
    path("*_reheader.bam.header")
    script:
    """
    bash reheader_chr.sh ${bamfile} ${meta.sample}_reheader.bam
    samtools flagstat ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.flagstat
    samtools view -H ${meta.sample}_reheader.bam > ${meta.sample}_reheader.bam.header
    """
}

    //samtools view -H ${bamfile} | sed '/^@SQ/s/SN\:chr/SN\:/' | samtools reheader - ${bamfile} > ${meta}_reheader.bam
    //samtools index ${meta}_reheader.bam
   // samtools flagstat ${meta}_reheader.bam > ${meta}_reheader.bam.flagstat
   // samtools view -H ${meta}_reheader.bam > ${meta}_reheader.bam.header

workflow alt_align{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes
    reference_basename
    chr
    main:
    bwa_mem_align_alt(
        ch_ref,
        ch_fastq,
        reference_basename
    )
    bwa_mem_postalt(
        ch_ref,
        bwa_mem_align_alt.out.samfile,
        reference_basename
    )
    /*
    bwa_mem_align_alt_postalt(
       ch_ref,
       ch_fastq,
       reference_basename
    ) 
    */  
    samtools_sort(
        bwa_mem_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    markduplicates(
        samtools_sorted_index
    )
    extractContigs(
        ch_hlatypes,
        ch_ref,
        reference_basename,
        chr
    )
    fasta_index_bed(
        ch_ref,
        chr,
        reference_basename
    )
    subsetBam(
        markduplicates.out.markdupbam,
        extractContigs.out.alt_contigs,
        extractContigs.out.hla_contigs,
        ch_ref,
        fasta_index_bed.out.fasta_bed
    )
    emit: subsetBam.out.subsetbam
}

workflow alt_align_chr19{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes
    reference_basename

    main:
    /*
    bwa_mem_align_alt(
        ch_ref,
        ch_fastq,
        reference_basename
    )
    bwa_mem_postalt(
        ch_ref,
        bwa_mem_align_alt.out.samfile,
        reference_basename
    )
    */
    bwa_mem_align_alt_postalt(
       ch_ref,
       ch_fastq,
       reference_basename
    )
    samtools_sort(
        bwa_mem_align_alt_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    markduplicates(
        samtools_sorted_index
    )
    extractContigsTest(
        ch_hlatypes,
        ch_ref,
        reference_basename
    )
    fasta_index_bed(
        ch_ref,
        "chr19",
        reference_basename
    )
    subsetBamtest(
        markduplicates.out.markdupbam,
        extractContigsTest.out.alt_chr19_contigs,
        extractContigsTest.out.hla_contigs,
        ch_ref,
        fasta_index_bed.out.fasta_bed
    )
    emit:
    subsetBamtest.out.subsetbam
}

workflow prepPolysolver{
    /*
    need to realign without alt contigs and remove "chr" from bam header
    */
    take: 
    subsetbam
    reference
    reference_basename

    main:
    bam2fastq(
        subsetbam
    )
    realignwithoutAlt(
        bam2fastq.out.subsetfastq,
        reference,
        reference_basename
        ) 
    samtools_sort(
        realignwithoutAlt.out.realignbam
        )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    samtools_sorted_index = samtools_sort.out.sortedbam.join(samtools_index.out.bam_indexed, by: 0)
    reheaderChr(
        samtools_sorted_index
    )
}

/*
workflow {
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
    //ch_fastq.view()
    reference_basename = Channel.value(params.reference_basename)
    //Channel
    //.fromPath(params.reference_dir, checkIfExists: true)
    //.set { ch_ref }
    ch_ref = file(params.reference_dir, checkIfExists: true)
    //ch_ref.view()
    //Channel
    //.fromPath(params.hlatypes, checkIfExists: true)
    //.set { ch_hlatypes }
    ///.view()
    ch_hlatypes = file(params.hlatypes, checkIfExists: true)
    chromosome = Channel.value(params.chr)
    alt_align(  
    ch_fastq,
    ch_ref,
    ch_hlatypes,
    reference_basename,
    chromosome
    )
    //alt_align_chr19.out.collect().view()
    prepPolysolver(
    alt_align.out.collect(),
    ch_ref,
    reference_basename
    )
}
*/

workflow {
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample')
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
    reference_basename = Channel.value(params.reference_basename)
    ch_ref = file(params.reference_dir, checkIfExists: true)
    ch_hlatypes = file(params.hlatypes, checkIfExists: true)
    chromosome = Channel.value(params.chr)
    alt_align(
    ch_fastq,
    ch_ref,
    ch_hlatypes,
    reference_basename,
    chromosome
    )
    prepPolysolver(
    alt_align.out,
    ch_ref,
    reference_basename
    )
}

