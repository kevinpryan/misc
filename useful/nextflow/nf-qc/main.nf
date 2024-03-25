process fastqc{
    publishDir "${params.outdir}/${meta.sample}/fastqc"

 input:
 tuple val(meta), path(reads)
 //val parstr

  output:
  tuple val(meta), path("*"), emit: fastqc_reports
  tuple val(meta), path("*zip"), emit: fastqc_zips
  tuple val(meta) path("*html"), emit: fastqc_htmls

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
    ch_fastq.view()
    //fastqc(ch_fastq)
    //multiqc(fastqc.fastqc_reports.collect())
}