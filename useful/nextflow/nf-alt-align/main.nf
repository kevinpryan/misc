#!/usr/bin/env nextflow
params.reads_1 = '/data2/caf_wes_bc/raw_data/*_{1,2}.fastq.gz'
params.reads_2 = '/data2/caf_wes_bc/raw_data_final_samples*_{1,2}.fastq.gz'

Channel
    .fromPath(params.reads_1)
    .ifEmpty { error "No input files found: $params.reads_1" }
    .set { input_ch }

Channel
    .fromPath(params.reads_2)
    .ifEmpty { error "No input files found: $params.reads_2" }
    .set { input_ch2 }

// merge the two channels
input_ch.merge(input_ch2)
    .set { input_ch_merged }

input_ch_merged.view() 