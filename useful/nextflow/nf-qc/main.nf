process fastqc{
    publishDir "${params.outdir}/${meta.sample}/fastqc"

 input:
 tuple val(meta), path(reads)
 //val parstr

  output:
  tuple val(meta), path("*"), emit: fastqc_reports
  tuple val(meta), path("*zip"), emit: fastqc_zips
  tuple val(meta), path("*html"), emit: fastqc_htmls
  path("*zip"), emit: fastqc_zips_path
  script:
  """
  fastqc ${reads} -t ${task.cpus}
  """
}

process multiqc{
    publishDir "${params.outdir}/multiqc"
    input:
    path metric_files
   
    output:
    path "multiqc*", emit: multiqc_reports

    script:
    """
    multiqc .
    """
}

process strandedness{
    publishDir "${params.outdir}/${meta.sample}/strandedness"
    input:
    tuple val(meta), path(reads)
    path gtf
    output:
    path "*", emit: strandedness
    script:
    """
    check_strandedness -g ${gtf} -fa /data/kryan/reference/gencode.v40.transcripts.fa -r1 $read1 -r2 $read2 -k /data/kryan/rna_seq_bc/caf_subtypes/EGAD00001005744/kallisto_index > ${sample}_strandedness.txt
    strandedness ${reads}
    """

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
        ch_ref = file(params.reference_dir, checkIfExists: true)
    fastqc(ch_fastq)
    multiqc(fastqc.out.fastqc_zips_path.collect())
}
