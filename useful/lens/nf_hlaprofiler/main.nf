#!/usr/bin/env nextflow

process get_fastqs {
// Symlink and emit FASTQ pairs for a FASTQ Prefix of a Patient Name from a Dataset
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(prefix) - FASTQ Prefix
//     val(dataset) - Dataset
//     val(run_name) - Run Name
//     val fq_dir - Project-specific /fastq directory
//
// output:
//   tuple => emit: fastqs
//     val(pat_name) - Patient Name
//     val(run_name) - Sequence run name. Unique within dataset.
//     val(dataset) - Dataset
//     path("${prefix}*1*.f*q.gz") - FASTQ 1
//     path("${prefix}*2*.f*q.gz") - FASTQ 2

  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(run_name), val(dataset) //, val(prefix)
  val fq1
  val fq2
  //val fq_dir

  output:
  tuple val(pat_name), val(run_name), val(dataset), path("*${run_name}_1*.f*q.gz"), path("*${run_name}_2*.f*q.gz"), emit: fastqs

  script:
  """
  fastq1=${fq1}
  fastq2=${fq2}
  """
}

//fastq1=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_1\\|R1\\).*\\.\\(fq\\|fastq\\)\\.gz")
//fastq2=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_2\\|R2\\).*\\.\\(fq\\|fastq\\)\\.gz")
//ln -s \${fastq1} ${dataset}-${pat_name}-${run_name}_1.fastq.gz
//ln -s \${fastq2} ${dataset}-${pat_name}-${run_name}_2.fastq.gz

process seqtk_sample {
// require:
//   FQS
//   params.seqtk$seqtk_sample_read_count
//   SEED
//   params.seqtk$seqtk_sample_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'seqtk_container'
  label 'seqtk_sample'
  publishDir "${params.outdir}/${dataset}/${pat_name}/${run}/seqtk_sample"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val count
  val seed
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*_1.*.subd.*.fastq.gz"), path("*_2.*.subd.*.fastq.gz"), emit: subd_fqs

  script:
  """
  FQ1_BUFFER=`echo ${fq1}`
  FQ2_BUFFER=`echo ${fq2}`

  FQ1_OUT=`echo \${FQ1_BUFFER} | sed 's/.fastq.gz//g' | sed 's/.fq.gz//g'`
  FQ2_OUT=`echo \${FQ2_BUFFER} | sed 's/.fastq.gz//g' | sed 's/.fq.gz//g'`

  HR=`numfmt --to=si ${count}`

  seqtk sample -s${seed} ${parstr} ${fq1} ${count} | gzip -c > \${FQ1_OUT}${suffix}.subd.\${HR}.fastq.gz&
  seqtk sample -s${seed} ${parstr} ${fq2} ${count} | gzip -c > \${FQ2_OUT}${suffix}.subd.\${HR}.fastq.gz
  """
}

process trim_galore_hlap {
// Runs trim galore with hard trimming for HLAProfiler. Needs to be wrapped up
// into the primary trim_galore process definition.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - FASTQ Prefix
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   val parstr - Parameter String
//
// output:
//   tuple => emit: procd_fqs
//     val(pat_name) - Patient Name
//     val(run) - FASTQ Prefix
//     val(dataset) - Dataset
//     path("${dataset}-${pat_name}-${run}*_1*.trimmed.f*q.gz") - Trimmed FASTQ 1
//     path("${dataset}-${pat_name}-${run}*_2*.trimmed.f*q.gz") - Trimmed FASTQ 2
//   path("meta") - Metadata File

// require:
//   FQS
//   params.trim_galore$trim_galore_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'trim_galore_container'
  label 'trim_galore'
//  publishDir "${params.samps_out_dir}/${dataset}/${run}/trim_galore"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*_1*.trimmed*50bp_5prime.f*q.gz"), path("${dataset}-${pat_name}-${run}*_2*.trimmed*50bp_5prime.f*q.gz"), emit: procd_fqs

  script:
  """
  trim_galore --basename ${run} ${parstr} ${fq1} ${fq2} -j ${task.cpus}
  """
}

process hlaprofiler_predict {
// Runs HLAProfiler.pl predict
// Uses the database included with the Docker image.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: calls
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*HLATypes.txt") - Calls file
//
// require:
//   FQS
//   params.hlaprofiler$hlaprofiler_predict$parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'hlaprofiler_container'
  label 'hlaprofiler_predict'
  publishDir "${params.outdir}/${dataset}/${pat_name}/${run}/hlaprofiler_predict"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  //val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*HLATypes.txt"), emit: calls

  script:
  """
  perl /HLAProfiler-v.1.10beta/bin/HLAProfiler.pl predict \
    -fastq1 ${fq1} \
    -fastq2 ${fq2} \
    -database_name hla_database \
    -database_dir /opt/HLAProfiler-1.0.0-db_only \
    -reference /opt/HLAProfiler-1.0.0-db_only/hla_database/data/reference/hla.ref.merged.fa \
    -kraken_path /HLAProfiler-v.1.10beta/kraken-0.10.5-beta-ea.1 \
    -output_dir subwork \
    -allele_refinement all \
    -l sample.HLAProfiler.log \
    -threads ${task.cpus}
  mv subwork/*/*HLATypes.txt ${dataset}-${pat_name}-${run}.HLATypes.txt
  find subwork -name "*fq" -exec gzip {} \\;
  """
}

//     ${parstr} \
dataset = Channel.of(params.dataset)
run = Channel.of(params.run)
pat_name = Channel.of(params.pat_name)
ch_fq1 = Channel.fromPath(params.fq1, checkIfExists: true)
ch_fq2 = Channel.fromPath(params.fq2, checkIfExists: true)
seqtk_read_count = Channel.of(params.seqtk_read_count)
seqtk_seed = Channel.of(params.seqtk_seed)
seqtk_sample_suffix = Channel.of(params.seqtk_sample_suffix)
trim_galore_hlap_parameters = Channel.of(params.trim_galore_hlap_parameters)
hlaprofiler_parameters = Channel.of(params.hlaprofiler_parameter)

/*
workflow procd_fqs_to_hlaprofiler_calls {
// require:
//   FQS
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_read_count
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_seed
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_sample_suffix
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_parameters
//   params.immuno$procd_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters
//   params.immuno$procd_fqs_to_hlaprofiler_calls$hlaprofiler_parameters
  take:
    fqs
    seqtk_read_count
    seqtk_seed
    seqtk_sample_suffix
    seqtk_parameters
    trim_galore_hlap_parameters
    hlaprofiler_parameters
  main:
    seqtk_sample(
      fqs,
      seqtk_read_count,
      seqtk_seed,
      seqtk_sample_suffix,
      seqtk_parameters)
    trim_galore_hlap(
      seqtk_sample.out.subd_fqs,
      trim_galore_hlap_parameters) 
    hlaprofiler_predict(
      trim_galore_hlap.out.procd_fqs,
      hlaprofiler_parameters)
  emit:
    calls = hlaprofiler_predict.out.calls
}
*/

workflow {
    get_fastqs_input = pat_name.concat(run_name, dataset)
    get_fastqs_input.view()
    )
}