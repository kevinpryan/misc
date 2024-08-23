#!/usr/bin/env nextflow

process fastqc{
    publishDir "${params.outdir}/${meta.sample}/fastqc"
  conda (params.enable_conda ? "bioconda::fastqc=0.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    }
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
conda (params.enable_conda ? "bioconda::multiqc=1.24" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.24--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.24--pyhdfd78af_0"
    }
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
   conda (params.enable_conda ? "bioconda::how_are_we_stranded_here=1.0.1" : null)
   if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/how_are_we_stranded_here:1.0.1--pyhfa5458b_0"
    } else {
        container "quay.io/biocontainers/how_are_we_stranded_here:1.0.1--pyhfa5458b_0"
    }
    publishDir "${params.outdir}/strandedness"
    input:
    tuple val(meta), path(reads)
    path gtf
    path fasta
    output:
    path "*_strandedness.txt", emit: strandedness
    script:
    """
    check_strandedness -g ${gtf} -fa ${fasta} -r1 *1.fastq.gz -r2 *2.fastq.gz > ${meta.sample}_strandedness.txt
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
    //Channel.fromPath(params.samplesheet, checkIfExists: true)
    ///| splitCsv( header:true )
    //| map { row ->
    //    meta = row.subMap('sample')
    //    [meta, [
    ///        file(row.fastq_1, checkIfExists: true),
    //        file(row.fastq_2, checkIfExists: true)]]
    //}
   // | groupTuple(by: 1)
    //| set { ch_fastq_tuple }
    //ch_fastq_tuple.view()

    ch_gtf = file(params.gtf, checkIfExists: true)
    ch_fasta = file(params.fasta, checkIfExists: true)
    strandedness(ch_fastq, 
                ch_gtf,
                ch_fasta)
    fastqc(ch_fastq)
    multiqc(fastqc.out.fastqc_zips_path.collect())
}
