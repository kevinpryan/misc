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