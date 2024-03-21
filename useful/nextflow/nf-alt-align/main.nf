#!/usr/bin/env nextflow

process bwa_mem_align_alt{
    publishDir "$params.outdir/align"

    input:
    path reference
    tuple val(meta), path(fq1), path(fq2)

    output:
    tuple val(meta), path("*.bwamem.sam"), emit: samfile
    path("*.bwamem.sam.header")
    path("*.bwamem.sam.flagstat")
    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${fq1} ${fq2} > "${meta}.bwamem.sam"
    samtools view -H ${meta}.bwamem.sam > ${meta}.bwamem.sam.header
    samtools flagstat ${meta}.bwamem.sam > ${meta}.bwamem.sam.flagstat
    """
}

process bwa_mem_postalt{
    publishDir "$params.outdir/postalt"

    input:
    path reference
    tuple val(meta), path(samfile)

    output:
    tuple val(meta), path("*_postalt.bam"), emit: bamfile_postalt
    path("*_postalt.sam.header")
    path("*_postalt.sam.flagstat")
    path("*_postalt.bam.flagstat")

    script:
    """
    k8-linux bwa-postalt.js \
    hs38DH.fa.alt \
    ${meta}.bwamem.sam > ${meta}_postalt.sam
    samtools view -H ${meta}_postalt.sam > ${meta}_postalt.sam.header
    samtools flagstat ${meta}_postalt.sam > ${meta}_postalt.sam.flagstat
    samtools view -bh -o ${meta}_postalt.bam ${meta}_postalt.sam
    samtools flagstat ${meta}_postalt.bam > ${meta}_postalt.bam.flagstat
    """
}

process samtools_sort{
    publishDir "$params.outdir/sort"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sortedbam

    script:
    """
    samtools sort -o ${meta}.sorted.bam ${bamfile}
    """
}

process samtools_index{
    publishDir "$params.outdir/index"

    input:
    tuple val(meta), path(sortedbam)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam_indexed

    script:
    """
    samtools index ${sortedbam}
    """
}   

process markduplicates{
    publishDir "$params.outdir/markduplicates"

    input:
    tuple val(meta), path(sortedbam)

    output:
    tuple val(meta), path("*_sorted_mdup.bam"), emit: markdupbam
    path("*.sorted.markdup.bam.flagstat")

    script:
    """
    bammarkduplicates I=${sortedbam} O=${meta}_sorted_mdup.bam" index=1 rmdup=0
    samtools flagstat ${meta}_sorted_mdup.bam > ${meta}_sorted_mdup.bam.flagstat
    """
}

process extractContigs {
    publishDir "$params.outdir/extractContigs"

    input:
    path hlatypes
    path reference
    output:
    path("bwakit-alt_chr6_contigs.txt"), emit: alt_chr6_contigs
    path("bwakit-alt-chr19_contigs.txt"), emit: alt_chr19_contigs
    path("bwakit-hla_contigs.txt"), emit: hla_contigs
    script:
    """
    grep 'chr6_' hs38DH.fa.alt | awk -F '\\t' '{print \$1}' > bwakit-alt_chr6_contigs.txt
    awk -F '\\t' '{print \$3}' ${hlatypes} > bwakit-hla_contigs.txt
    """
}

process subsetBam{
    publishDir "$params.outdir/subsetBam"

    input:
    tuple val(meta), path(bamfile)
    path alt_chr6_contigs
    path hla_contigs
    output:
    tuple val(meta), path("*_subset.sorted.bam"), path("*_subset.sorted.bam.bai"), path("*_subset.sorted.bam.flagstat"), path("*_subset.sorted.bam.flagstat"), emit: subsetbam

    script:
    """
    primary_contig='chr6:28509970-33480727'
    no_contig='*'
    samtools view -o ${meta}_subset.bam -b ${bamfile} $primary_contig ${alt_chr6_contigs} ${hla_contigs} "$no_contig" 
    samtools sort -o ${meta}_subset.sorted.bam ${meta}_subset.bam
    sambamba index ${meta}_subset.sorted.bam
    samtools view -H ${meta}_subset.sorted.bam > ${meta}_subset.sorted.bam.header
    samtools flagstat ${meta}_subset.sorted.bam > ${meta}_subset.sorted.bam.flagstat
    """
}

process subsetBamtest{
    publishDir "$params.outdir/subsetBam"

    input:
    tuple val(meta), path(bamfile)
    path alt_chr19_contigs
    path hla_contigs
    output:
    tuple val(meta), path("*_subset.sorted.bam"), path("*_subset.sorted.bam.bai"), path("*_subset.sorted.bam.flagstat"), path("*_subset.sorted.bam.flagstat"), emit: subsetbam

    script:
    """
    primary_contig='chr6:28509970-33480727'
    no_contig='*'
    samtools view -o ${meta}_subset.bam -b ${bamfile} $primary_contig ${alt_chr19_contigs} ${hla_contigs} "$no_contig" 
    samtools sort -o ${meta}_subset.sorted.bam ${meta}_subset.bam
    sambamba index ${meta}_subset.sorted.bam
    samtools view -H ${meta}_subset.sorted.bam > ${meta}_subset.sorted.bam.header
    samtools flagstat ${meta}_subset.sorted.bam > ${meta}_subset.sorted.bam.flagstat
    """
}

process bam2fastq{
    publishDir "$params.outdir/bam2fastq"

    input:
    tuple val(meta), path(subsetbam)
    output:
    tuple val(meta), path("*_subset.1.fq"), path("*_subset.2.fq"), emit: subsetfastq

    script:
    """
    sambamba view \
    -f "bam" -h -p -l 0 -t ${task.cpus} \
    "$IN" | sambamba sort -p -n -t ${task.cpus} -o - /dev/stdin |
    samtools fastq -@ 32 /dev/stdin \
    -1 "${meta}_subset.1.fq" \
    -2 "${meta}_subset.2.fq" \
    -0 /dev/null -s /dev/null -n
    """
}

process realignwithoutAlt{
    publishDir "$params.outdir/realignwithoutAlt"
    
    input:
    tuple val(meta), path(fq1), path(fq2)
    path reference

    output:
    tuple val(meta), path("*_realign.sorted.bam"), path("*_realign.sorted.bam.bai"), path("*_realign.sorted.bam.flagstat"), path("*_realign.sorted.bam.flagstat"), emit: realignbam

    script:
    """
    bwa mem -t ${task.cpus} -j ${reference} ${fq1} ${fq2} > ${meta}_realign.sam
    samtools view -H ${meta}_realign.sam > ${meta}_realign.sam.header
    samtools flagstat ${meta}_realign.sam > ${meta}_realign.sam.flagstat
    """
}

process reheaderChr{
    publishDir "$params.outdir/reheaderChr"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*_reheader.bam"), path("*_reheader.bam.bai"), path("*_reheader.bam.flagstat"), path("*_reheader.bam.flagstat"), emit: reheaderbam

    script:
    """
    sed 's/SN:chr/SN:/' <(${samtools.view("-H", bamfile)}) | ${samtools.reheader("-", bamfile)} > ${meta}_reheader.bam
    ${samtools.index(meta + "_reheader.bam")}
    ${samtools.flagstat(meta + "_reheader.bam")} > ${meta}_reheader.bam.flagstat
    ${samtools.view("-H", meta + "_reheader.bam")} > ${meta}_reheader.bam.header
    """
}

    //samtools view -H ${bamfile} | sed '/^@SQ/s/SN\:chr/SN\:/' | samtools reheader - ${bamfile} > ${meta}_reheader.bam
    //samtools index ${meta}_reheader.bam
   // samtools flagstat ${meta}_reheader.bam > ${meta}_reheader.bam.flagstat
   // samtools view -H ${meta}_reheader.bam > ${meta}_reheader.bam.header

workflow alt_align_chr6{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes

    main:
    bwa_mem_align_alt(
        ch_ref,
        ch_fastq
    )
    bwa_mem_postalt(
        ch_ref,
        bwa_mem_align_alt.out.samfile
    )
    samtools_sort(
        bwa_mem_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    markduplicates(
        samtools_index.out.bam_indexed
    )
    extractContigs(
        ch_hlatypes,
        ch_ref
    )
    subsetBam(
        markduplicates.out.markdupbam,
        extractContigs.out.alt_chr6_contigs,
        extractContigs.out.hla_contigs
    )
}

workflow alt_align_chr19{
    take: 
    ch_fastq
    ch_ref
    ch_hlatypes

    main:
    bwa_mem_align_alt(
        ch_ref,
        ch_fastq
    )
    bwa_mem_postalt(
        ch_ref,
        bwa_mem_align_alt.out.samfile
    )
    samtools_sort(
        bwa_mem_postalt.out.bamfile_postalt
    )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    markduplicates(
        samtools_index.out.bam_indexed
    )
    extractContigs(
        ch_hlatypes,
        ch_ref
    )
    subsetBamtest(
        markduplicates.out.markdupbam,
        extractContigs.out.alt_chr6_contigs,
        extractContigs.out.hla_contigs
    )
}

workflow prepPolysolver{
    /*
    need to realign without alt contigs and remove "chr" from bam header
    */
    take: 
    subsetbam

    main:
    bam2fastq(
        subsetbam
    )
    realignwithoutAlt(
        bam2fastq.out.subsetfastq,
        reference
        ) 
    samtools_sort(
        relalignwithoutAlt.out.realignbam
        )
    samtools_index(
        samtools_sort.out.sortedbam
    )
    reheaderChr(
        samtools_index.out.bam_indexed
    )
}


workflow {
    Channel.fromPath(params.samplesheet)
    | splitCsv( header:true )
    | map { row ->
        meta = [id:row.sample]
        [meta, [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)]]
    }
    | set { ch_fastq }
    Channel
    .fromPath(params.reference_dir)
    //.ifEmpty { error "No reference files found: $params.reference_dir" }
    .set { ch_ref }
    .view()
    Channel
    .fromPath(params.hlatypes)
    //.ifEmpty { error "No hla types found: $params.hlatypes" }
    .set { ch_hlatypes }
    .view()
    //alt_align_chr19(ch_fastq, ch_ref)

}