#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { VCF2BED } from './modules/reusable'
include { SUBSET_CRAM as SUBSET_CRAM_CAF } from './modules/reusable'
include { SUBSET_CRAM as SUBSET_CRAM_TAN } from './modules/reusable'

// Read in the list of VCF files

vcf_files_caf = Channel.fromPath(params.vcf_path)
    .map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[0], it, it.simpleName ] }

vcf_files_tan = Channel.fromPath(params.vcf_path)
    .map { it -> [ it.simpleName.toString().split('\\.')[0].split('\\_')[2], it, it.simpleName ] }
// Read in the list of CRAM files
cram_files = Channel.fromPath(params.cram_path)
     .map{ it -> [ it.simpleName.toString(), it]}
crai_files = Channel.fromPath(params.crai_path)

cram_crai = cram_files
            .join(crai_files)
// Combine the VCF and CRAM files based on the sample IDs
vcf_cram_pairs_caf = cram_files
                     .join(vcf_files_caf)

vcf_cram_pairs_tan = cram_files
                     .join(vcf_files_tan)

vcf_cram_combined_full = vcf_cram_pairs_caf
                         .join(vcf_cram_pairs_tan, by:3)
                         .unique()

norm_vcfs_input = vcf_cram_pairs_caf
                  .map { it -> [ it[3], it[2]]}
                  .groupTuple()

ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
                  .collect()

workflow {
    VCF2BED(norm_vcfs_input)
    ch_bed_cram = VCF2BED.out.bed_from_vcf.join(vcf_cram_combined_full)
    ch_bed_cram_caf = ch_bed_cram.map { it -> [it[0], it[3], it[1]] }.groupTuple()
    ch_bed_cram_tan = ch_bed_cram.map { it -> [ it[0], it[6], it[1]]}.groupTuple()
    SUBSET_CRAM_CAF(ch_bed_cram_caf, ch_fasta)
    SUBSET_CRAM_TAN(ch_bed_cram_tan, ch_fasta)
}
