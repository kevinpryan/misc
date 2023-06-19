include { optitype } from 'modules/optitype/optitype.nf'

params.fq_dir = "${params.project_dir}/fastqs"
params.manifest = "${params.project_dir}/files/manifest.tsv"

ch_fq_dir = Channel.fromPath(params.fq_dir, checkIfExists: true)
ch_manifest = Channel.fromPath(params.manifest, checkIfExists: true)

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
  tuple val(pat_name), val(run_name), val(dataset), val(prefix)
  val fq_dir

  output:
  tuple val(pat_name), val(run_name), val(dataset), path("*${run_name}_1*.f*q.gz"), path("*${run_name}_2*.f*q.gz"), emit: fastqs

  script:
  """
  fastq1=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_1\\|R1\\).*\\.\\(fq\\|fastq\\)\\.gz")
  fastq2=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_2\\|R2\\).*\\.\\(fq\\|fastq\\)\\.gz")
  ln -s \${fastq1} ${dataset}-${pat_name}-${run_name}_1.fastq.gz
  ln -s \${fastq2} ${dataset}-${pat_name}-${run_name}_2.fastq.gz
  """
}

workflow manifest_to_raw_fqs {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest.map{ [it[0], it[1], it[2], it[3]] }, //tuple val(pat_name), val(prefix), val(dataset), val(run)
      params.fq_dir)
  emit:
    fastqs = get_fastqs.out.fastqs
}

// not actually run, don't need to trim to use optitype
workflow manifest_to_rna_procd_fqs {
// require:
//   MANIFEST
//   params.manifest_to_rna_procd_fqs$trim_galore_parameters
  take:
    manifest
    //trim_galore_parameters
  main:
    manifest_to_raw_fqs(
      manifest.filter{ it[4] =~ /RNA/ })
    trim_galore(
      manifest_to_raw_fqs.out.fastqs,
      trim_galore_parameters
    )
  emit:
    raw_fqs = manifest_to_raw_fqs.out.fastqs
    procd_fqs = trim_galore.out.procd_fqs
    fastqc_zips = trim_galore.out.fastqc_zips
}


workflow procd_fqs_to_optitype_calls {
// require:
//   FQS
  take:
    fqs
  main:
    optitype(
      fqs)
  emit:
    calls = optitype.out.calls
}

workflow {
    manifest_to_raw_fqs(
        ch_manifest
    )
    procd_fqs_to_optitype_calls(
        manifest_to_raw_fqs.out.fastqs
    )
}
