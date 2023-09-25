#!/usr/bin/env nextflow

// Preprocessing
//include { manifest_to_alns as manifest_to_dna_alns } from '../alignment/alignment.nf'
include { manifest_to_alns as manifest_to_rna_alns } from '../alignment/alignment.nf'

// Sorting and indexing transcriptome alignment
include { samtools_sort as samtools_sort_txome } from  '../samtools/samtools.nf'
include { samtools_index as samtools_index_txome } from  '../samtools/samtools.nf'

// Indexing RNA (genome) alignment
include { samtools_index as samtools_index_rna } from  '../samtools/samtools.nf'

// Indexing viral alignment
//include { samtools_index as samtools_index_viral } from  '../samtools/samtools.nf'

// BAM sanitization
//include { alns_to_procd_alns as alns_to_dna_procd_alns } from '../alignment/alignment.nf'
include { alns_to_procd_alns as alns_to_rna_procd_alns } from '../alignment/alignment.nf'

// Indexing DNA alignment
//include { samtools_index as samtools_index_dna } from  '../samtools/samtools.nf'

// Transcript quantificaiton
include { alns_to_transcript_counts } from '../rna_quant/rna_quant.nf'

// Somatic variant calling, filtering, intersecting, and annotating
//include { alns_to_som_vars } from '../somatic/somatic.nf'/
//include { som_vars_to_filtd_som_vars } from '../somatic/somatic.nf'
//include { som_vars_to_normd_som_vars } from '../somatic/somatic.nf'
//include { som_vars_to_isecd_som_vars } from '../somatic/somatic.nf'
//include { som_vars_to_union_som_vars } from '../somatic/somatic.nf'
//include { htslib_bgzip_somatic } from '../htslib/htslib.nf'
//include { htslib_bgzip_somatic as htslib_bgzip_somatic_isec } from '../htslib/htslib.nf'
//include { snpeff_ann } from '../snpeff/snpeff.nf'
//i//nclude { snpsift_filter as snpsift_filter_snvs } from '../snpeff/snpeff.nf'
//include { snpsift_filter as snpsift_filter_indels} from '../snpeff/snpeff.nf'

// Germline variant calling and filtering
//include { alns_to_germ_vars } from '../germline/germline.nf'
//include { germ_vars_to_filtd_germ_vars } from '../germline/germline.nf'
//include { htslib_bgzip } from '../htslib/htslib.nf'

// Creating tumor VCF (germline + somatic variants)
//include { bcftools_index } from '../bcftools/bcftools.nf'
//i//nclude { bcftools_index_somatic } from '../bcftools/bcftools.nf'
//include { germ_and_som_vars_to_tumor_vars } from '../seq_variation/seq_variation.nf'

// Phasing of variants
//include { make_phased_tumor_vars } from '../seq_variation/seq_variation.nf'
//include { make_phased_germline_vars } from '../seq_variation/seq_variation.nf'

// Viral expression detection
//include { alns_to_viruses } from '../viral/viral.nf'
//include { unaligned_fqs_to_virdetect_cds_counts } from '../viral/viral.nf'

// Splice variant detection
//include { alns_to_splice_variants } from '../splice/splice.nf'

// Fusion detection
//include { procd_fqs_to_fusions } from '../fusion/fusion.nf'
/*
// Neoantigen peptides generation
include { som_vars_to_neos as snvs_to_neos} from '../neos/neos.nf'
include { som_vars_to_neos as indels_to_neos} from '../neos/neos.nf'
include { ervs_to_neos } from '../neos/neos.nf'
include { selfs_to_neos } from '../neos/neos.nf'
include { fusions_to_neos } from '../neos/neos.nf'
include { viruses_to_neos } from '../neos/neos.nf'
*/
/*
// Join peptide fastas
include { combine_peptide_fastas } from '../immuno/immuno.nf'
include { combine_nt_fastas } from '../immuno/immuno.nf'

// MHC calling
include { user_provided_alleles_to_netmhcpan_alleles } from '../immuno/immuno.nf'
include { extract_alleles_from_manifest } from '../immuno/immuno.nf'
include { procd_fqs_to_mhc_alleles } from '../immuno/immuno.nf'

// TCR repertoire
include { procd_fqs_to_tcr_repertoire } from '../immuno/immuno.nf'

// CNAs/CCF
include { alns_to_cnas } from '../onco/onco.nf'
include { cnas_and_vcfs_to_ccfs } from '../onco/onco.nf'

// pMHC summarization
include { peps_and_alleles_to_antigen_stats } from '../immuno/immuno.nf'

// Aggregate pMHC summaries
include { aggregate_pmhc_summaries } from '../immuno/immuno.nf'

// Annotate pMHC summaries
include { lenstools_annotate_pmhcs } from '../lenstools/lenstools.nf'

// Agretopicity
include { calculate_agretopicity } from '../immuno/immuno.nf'

// Quantify pMHC CDS abunance
include { lenstools_get_erv_and_cta_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_viral_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_snv_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_indel_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_fusion_peptide_read_count } from '../lenstools/lenstools.nf'
// Needed read extraction process for fusion quantification
include { fusions_to_extracted_reads } from '../fusion/fusion.nf'
include { seqtk_subseq } from '../seqtk/seqtk.nf'
include { lenstools_get_splice_peptide_read_count } from '../lenstools/lenstools.nf'

include { lenstools_combine_read_counts } from '../lenstools/lenstools.nf'

// Report annotations
include { lenstools_add_generic_annotation as annotate_ccfs } from '../lenstools/lenstools.nf'

// Visualization
include { lenstools_make_lens_bed } from '../lenstools/lenstools.nf'
include { igv_snapshot_automator } from '../igv_snapshot_automator/igv_snapshot_automator.nf'
include { bam_subsetter } from '../samtools/samtools.nf'

// Make final report
include { lenstools_make_lens_report } from '../lenstools/lenstools.nf'
*/
/* Non-workflow related intermediate file deletion */
//include { clean_work_files as clean_nonsorted_txome_bams } from '../utilities/utilities.nf'


workflow manifest_to_lens_rna_aln_quant {
// require:
//   MANIFEST
  take:
    manifest

  main:

    // Initial DNA alignments
    /*
    manifest_to_dna_alns(
      manifest.filter{ it[4] =~ /dna|DNA|wes|WES|wxs|WXS|WGS|wgs/ },
      params.lens$alignment$manifest_to_dna_alns$fq_trim_tool,
      params.lens$alignment$manifest_to_dna_alns$fq_trim_tool_parameters,
      params.lens$alignment$manifest_to_dna_alns$aln_tool,
      params.lens$alignment$manifest_to_dna_alns$aln_tool_parameters,
      params.lens$alignment$manifest_to_dna_alns$aln_ref,
      params.lens$alignment$manifest_to_dna_alns$gtf,
      '')
*/
    // Sanitizing DNA alignments
/*
    alns_to_dna_procd_alns(
      manifest_to_dna_alns.out.alns,
      '',
      '',
      params.lens$alignment$alns_to_dna_procd_alns$aln_ref,
      params.lens$alignment$alns_to_dna_procd_alns$bed,
      params.lens$alignment$alns_to_dna_procd_alns$gtf,
      params.lens$alignment$alns_to_dna_procd_alns$dup_marker_tool,
      params.lens$alignment$alns_to_dna_procd_alns$dup_marker_tool_parameters,
      params.lens$alignment$alns_to_dna_procd_alns$base_recalibrator_tool,
      params.lens$alignment$alns_to_dna_procd_alns$base_recalibrator_tool_parameters,
      params.lens$alignment$alns_to_dna_procd_alns$indel_realign_tool,
      params.lens$alignment$alns_to_dna_procd_alns$indel_realign_tool_parameters,
      params.lens$alignment$alns_to_dna_procd_alns$known_sites_ref,
      manifest.filter{ it[4] =~ /dna|DNA|wes|WES|wxs|WXS|WGS|wgs/ })
*/
    // Initial RNA alignments
    manifest_to_rna_alns(
      manifest.filter{ it[4] =~ /RNA/ },
      params.lens$alignment$manifest_to_rna_alns$fq_trim_tool,
      params.lens$alignment$manifest_to_rna_alns$fq_trim_tool_parameters,
      params.lens$alignment$manifest_to_rna_alns$aln_tool,
      params.lens$alignment$manifest_to_rna_alns$aln_tool_parameters,
      params.lens$alignment$manifest_to_rna_alns$aln_ref,
      params.lens$alignment$manifest_to_rna_alns$gtf,
      '')

/*
    // Somatic variant calling
    alns_to_som_vars(
      alns_to_dna_procd_alns.out.procd_bams,
      params.lens$somatic$alns_to_som_vars$som_var_caller,
      params.lens$somatic$alns_to_som_vars$som_var_caller_parameters,
      params.lens$somatic$alns_to_som_vars$som_var_caller_suffix,
      params.lens$somatic$alns_to_som_vars$aln_ref,
      params.lens$somatic$alns_to_som_vars$bed,
      params.lens$somatic$alns_to_som_vars$som_var_pon_vcf,
      params.lens$somatic$alns_to_som_vars$som_var_af_vcf,
      params.lens$somatic$alns_to_som_vars$known_sites_ref,
      params.lens$somatic$alns_to_som_vars$species,
      manifest.filter{ it[4] =~ /WES|wes|DNA|dna|WXS|wxs/})

    // Somatic variant filtering
    som_vars_to_filtd_som_vars(
      alns_to_som_vars.out.som_vars,
      params.lens$somatic$som_vars_to_filtd_som_vars$vcf_filtering_tool,
      params.lens$somatic$som_vars_to_filtd_som_vars$vcf_filtering_tool_parameters)

    // Somatic variant indel normalizing
    som_vars_to_normd_som_vars(
      som_vars_to_filtd_som_vars.out.filtd_som_vars,
      params.lens$somatic$som_vars_to_normd_som_vars$vcf_norming_tool,
      params.lens$somatic$som_vars_to_normd_som_vars$vcf_norming_tool_parameters,
      params.lens$somatic$som_vars_to_normd_som_vars$aln_ref)

    htslib_bgzip_somatic(
      som_vars_to_normd_som_vars.out.normd_som_vars)

    joint_vars = Channel.empty()

    if (params.lens$somatic$combine_strategy =~ /intersect/) {
      // Somatic variant intersectioning
        som_vars_to_isecd_som_vars(
          htslib_bgzip_somatic.out.bgzip_files,
          params.lens$somatic$som_vars_to_isecd_som_vars$vcf_isecing_tool,
          params.lens$somatic$som_vars_to_isecd_som_vars$vcf_isecing_tool_parameters)
        som_vars_to_isecd_som_vars.out.isecd_som_vars
          .set{ joint_vars }
     }
     if (params.lens$somatic$combine_strategy =~ /merge|union/) {
        som_vars_to_union_som_vars(
          htslib_bgzip_somatic.out.bgzip_files,
          params.lens$somatic$som_vars_to_union_som_vars$vcf_merging_tool,
          params.lens$somatic$som_vars_to_union_som_vars$vcf_merging_tool_parameters)
        som_vars_to_union_som_vars.out.union_som_vars
          .set{ joint_vars }
    }


    // Annotate somatic variants
    // This should be in an annot workflow that allows different tools.
    snpeff_ann(
//      som_vars_to_isecd_som_vars.out.isecd_som_vars,
      joint_vars,
      params.lens$snpeff$annot_tool_ref)


    // Filter somatic SNVs
    snpsift_filter_snvs(
      snpeff_ann.out.annotd_vcfs,
      params.lens$snpsift_filter_snvs$snpsift_snv_filter_parameters,
      "sfilt.snvs")

    // Filter somatic InDels
    snpsift_filter_indels(
      snpeff_ann.out.annotd_vcfs,
      params.lens$snpsift_filter_indels$snpsift_indel_filter_parameters,
      "sfilt.indels")


    // Germline variant calling
    alns_to_germ_vars(
      alns_to_dna_procd_alns.out.procd_bams.filter{ it[1] =~ /nd-/ },
      params.lens$germline$alns_to_germ_vars$germ_var_caller,
      params.lens$germline$alns_to_germ_vars$germ_var_caller_parameters,
      params.lens$germline$alns_to_germ_vars$germ_var_caller_suffix,
      params.lens$germline$alns_to_germ_vars$aln_ref,
      params.lens$germline$alns_to_germ_vars$bed)

    // Germline variant filtering
    germ_vars_to_filtd_germ_vars(
      alns_to_germ_vars.out.germ_vars,
      params.lens$germline$germ_vars_to_filtd_germ_vars$vcf_filtering_tool,
      params.lens$germline$germ_vars_to_filtd_germ_vars$vcf_filtering_tool_parameters)

    // Combine somatic and germline variants and phase to create "tumor"
    // variants (targetable variants coupled with phased, neighboring germline
    // and somatic variants.)

    htslib_bgzip(
      germ_vars_to_filtd_germ_vars.out.filtd_germ_vars)
    bcftools_index(
      htslib_bgzip.out.bgzip_files,
      params.lens$bcftools$bcftools_index$bcftools_index_parameters)
    htslib_bgzip_somatic_isec(
      snpeff_ann.out.annotd_vcfs)
    bcftools_index_somatic(
      htslib_bgzip_somatic_isec.out.bgzip_files,
      params.lens$bcftools$bcftools_index_somatic$bcftools_index_somatic_parameters)
    germ_and_som_vars_to_tumor_vars(
      bcftools_index.out.vcfs_w_csis,
      bcftools_index_somatic.out.vcfs_w_csis,
      params.lens$seq_variation$germ_and_som_vars_to_tumor_vars$vcf_merge_tool,
      params.lens$seq_variation$germ_and_som_vars_to_tumor_vars$vcf_merge_tool_parameters)
*/
/*
    // Sanitizing RNA alignments
    alns_to_rna_procd_alns(
      manifest_to_rna_alns.out.alns.filter{ it[1] =~ /ar-/ },
      manifest_to_rna_alns.out.junctions.filter{ it[1] =~ /ar-/},
      germ_and_som_vars_to_tumor_vars.out.tumor_vars,
      params.lens$alignment$alns_to_rna_procd_alns$aln_ref,
      params.lens$alignment$alns_to_rna_procd_alns$bed,
      params.lens$alignment$alns_to_rna_procd_alns$gtf,
      params.lens$alignment$alns_to_rna_procd_alns$dup_marker_tool,
      params.lens$alignment$alns_to_rna_procd_alns$dup_marker_tool_parameters,
      params.lens$alignment$alns_to_rna_procd_alns$base_recalibrator_tool,
      params.lens$alignment$alns_to_rna_procd_alns$base_recalibrator_tool_parameters,
      params.lens$alignment$alns_to_rna_procd_alns$indel_realign_tool,
      params.lens$alignment$alns_to_rna_procd_alns$indel_realign_tool_parameters,
      params.lens$alignment$alns_to_rna_procd_alns$known_sites_ref,
      manifest.filter{ it[4] =~ /rna|RNA|RNA-Seq|RNA-seq/ })
*/
/*
    // Phasing of tumor variants
    make_phased_tumor_vars(
      germ_and_som_vars_to_tumor_vars.out.tumor_vars_and_idxs,
      alns_to_dna_procd_alns.out.procd_bams_and_bais,
      alns_to_rna_procd_alns.out.procd_bams_and_bais,
      params.lens$seq_variation$make_phased_tumor_vars$aln_ref,
      params.lens$seq_variation$make_phased_tumor_vars$gtf,
      params.lens$seq_variation$make_phased_tumor_vars$species,
      params.lens$seq_variation$make_phased_tumor_vars$var_phaser_tool,
      params.lens$seq_variation$make_phased_tumor_vars$var_phaser_tool_parameters)
*/

/*
    // Phasing of germline variants
    make_phased_germline_vars(
      bcftools_index.out.vcfs_w_csis,
      alns_to_dna_procd_alns.out.procd_bams_and_bais,
      params.lens$seq_variation$make_phased_germline_vars$aln_ref,
      params.lens$seq_variation$make_phased_germline_vars$gtf,
      params.lens$seq_variation$make_phased_germline_vars$var_phaser_tool,
      params.lens$seq_variation$make_phased_germline_vars$var_phaser_tool_parameters)
*/

    // Transcript counts
    alns_to_transcript_counts(
      manifest_to_rna_alns.out.alt_alns, // Transcriptome aligns do not have InDels, not realignmement needed.
      params.lens$rna_quant$alns_to_transcript_counts$rna_ref,
      params.lens$rna_quant$alns_to_transcript_counts$gtf,
      params.lens$rna_quant$alns_to_transcript_counts$tx_quant_tool,
      params.lens$rna_quant$alns_to_transcript_counts$tx_quant_tool_parameters)

    // Sorting and indexing txome BAM for downstream application
    samtools_sort_txome(
      manifest_to_rna_alns.out.alt_alns.filter{ it[1] =~ 'ar-' },
      '')
    samtools_index_txome(
      samtools_sort_txome.out.bams,
      '')

    // Indexing RNA BAM for downstream application
    samtools_index_rna(
      manifest_to_rna_alns.out.alns.filter{ it[1] =~ 'ar-' },
      '')

    // Cleaning non-sorted txome BAM intermediates
    manifest_to_rna_alns.out.alt_alns
      .concat(samtools_sort_txome.out.bams)
      .groupTuple(by: [0, 1, 2], size: 2)
      .flatten()
      .filter{ it =~ /toTranscriptome.out.bam$/ }
      .set { nonsorted_txome_bam_done_signal }
//    clean_nonsorted_txome_bams(
//      nonsorted_txome_bam_done_signal)


/*
    // MHC calling
    pats_missing_alleles = Channel.empty()
    pats_with_alleles = Channel.empty()
    pats_with_standard_alleles = Channel.empty()
    extract_alleles_from_manifest(
        params.mhc_manifest) // Make this somehow accessible within workflow
    extract_alleles_from_manifest.out.patient_alleles.filter{ it[3] =~ /NA|Null|null/ }.set{ pats_missing_alleles }
    extract_alleles_from_manifest.out.patient_alleles.filter{ !(it[3] =~ /NA|Null|null/) }.set{ pats_with_alleles }
    user_provided_alleles_to_netmhcpan_alleles(
        pats_with_alleles)
    user_provided_alleles_to_netmhcpan_alleles.out.netmhcpan_alleles
      .set{ pats_with_standard_alleles }
    procd_fqs_to_mhc_alleles(
      manifest_to_rna_alns.out.procd_fqs.filter{ it[1] =~ 'ar-' },
      params.lens$immuno$procd_fqs_to_mhc_alleles$aln_tool,
      params.lens$immuno$procd_fqs_to_mhc_alleles$aln_tool_parameters,
      params.lens$immuno$procd_fqs_to_mhc_alleles$aln_ref,
      params.lens$immuno$procd_fqs_to_mhc_alleles$mhc_caller_tool,
      params.lens$immuno$procd_fqs_to_mhc_alleles$mhc_caller_tool_parameters)
    procd_fqs_to_mhc_alleles.out.alleles
      .concat(pats_with_standard_alleles)
      .set{ all_pat_mhc_alleles }

    // TCR repertoire
    procd_fqs_to_tcr_repertoire(
      manifest_to_rna_alns.out.procd_fqs.filter{ it[1] =~ 'ar-' },
      params.lens$immuno$procd_fqs_to_tcr_repertoire$tcr_rep_tool,
      params.lens$immuno$procd_fqs_to_tcr_repertoire$tcr_rep_tool_paraneters)


    // Splice variant detection
    alns_to_splice_variants(
      manifest_to_rna_alns.out.alns,
      params.splice$alns_to_splice_variants$splice_var_caller,
      params.splice$alns_to_splice_variants$splice_var_caller_parameters,
      params.splice$alns_to_splice_variants$splice_var_caller_ref,
      params.splice$alns_to_splice_variants$aln_ref,
      all_pat_mhc_alleles.filter{ it[1] =~ 'ar-' },
//      procd_fqs_to_mhc_alleles.out.alleles.filter{ it[1] =~ 'ar-' },
      params.splice$alns_to_splice_variants$gtf,
      params.splice$alns_to_splice_variants$gff,
      params.splice$alns_to_splice_variants$species,
      manifest)

    // Tumor Viral expression
    alns_to_viruses(
      manifest_to_rna_alns.out.alns.filter{ it[1] =~ 'ar-' },
      params.lens$viral$alns_to_viruses$viral_workflow,
      params.lens$viral$alns_to_viruses$viral_workflow_parameters,
      params.lens$viral$alns_to_viruses$viral_ref)
    unaligned_fqs_to_virdetect_cds_counts(
      alns_to_viruses.out.unaligned_fqs,
      params.lens$viral$unaligned_fqs_to_virdetect_cds_counts$viral_cds_ref)


    // Fusion detection
    procd_fqs_to_fusions(
      params.lens$fusion$procd_fqs_to_fusions$fusion_tool,
      params.lens$fusion$procd_fqs_to_fusions$fusion_tool_parameters,
      params.lens$fusion$procd_fqs_to_fusions$fusion_ref,
      params.lens$fusion$procd_fqs_to_fusions$dna_ref,
      params.lens$fusion$procd_fqs_to_fusions$gtf,
      manifest_to_rna_alns.out.procd_fqs.filter{ it[1] =~ 'ar-' })

    // CNA and cancer cell fraction
    alns_to_cnas(
      alns_to_dna_procd_alns.out.procd_bams,
      params.lens$onco$alns_to_cnas$cna_tool,
      params.lens$onco$alns_to_cnas$cna_tool_parameters,
      params.lens$onco$alns_to_cnas$aln_ref,
      params.lens$onco$alns_to_cnas$bed,
      manifest)

    cnas_and_vcfs_to_ccfs(
      alns_to_cnas.out.cnas,
      params.lens$onco$alns_to_cnas$cna_tool, // Referred here since it was used to generate CNAs.
//      som_vars_to_isecd_som_vars.out.isecd_som_vars,
      joint_vars,
      alns_to_som_vars.out.mutect2_vars,
      params.lens$onco$cnas_and_vcfs_to_ccfs$ccf_tool,
      params.lens$onco$cnas_and_vcfs_to_ccfs$ccf_tool_parameters)


    // SNVs to peptides
    snvs_to_neos(
      snpsift_filter_snvs.out.filtd_vcfs,
      make_phased_tumor_vars.out.phased_vcfs,
      make_phased_germline_vars.out.phased_vcfs,
      alns_to_transcript_counts.out.quants,
      params.lens$neos$snvs_to_neos$gtf,
      params.lens$neos$snvs_to_neos$dna_ref,
      params.lens$neos$snvs_to_neos$pep_ref,
      params.lens$neos$snvs_to_neos$som_var_type,
      params.lens$neos$snvs_to_neos$lenstools_filter_expressed_variants_parameters,
      params.lens$neos$snvs_to_neos$bcftools_index_phased_germline_parameters,
      params.lens$neos$snvs_to_neos$bcftools_index_phased_tumor_parameters,
      params.lens$neos$snvs_to_neos$lenstools_get_expressed_transcripts_bed_parameters,
      params.lens$neos$snvs_to_neos$samtools_faidx_fetch_somatic_folder_parameters)

    // InDels to peptides
    indels_to_neos(
      snpsift_filter_indels.out.filtd_vcfs,
      make_phased_tumor_vars.out.phased_vcfs,
      make_phased_germline_vars.out.phased_vcfs,
      alns_to_transcript_counts.out.quants,
      params.lens$neos$indels_to_neos$gtf,
      params.lens$neos$indels_to_neos$dna_ref,
      params.lens$neos$indels_to_neos$pep_ref,
      params.lens$neos$indels_to_neos$som_var_type,
      params.lens$neos$indels_to_neos$lenstools_filter_expressed_variants_parameters,
      params.lens$neos$indels_to_neos$bcftools_index_phased_germline_parameters,
      params.lens$neos$indels_to_neos$bcftools_index_phased_tumor_parameters,
      params.lens$neos$indels_to_neos$lenstools_get_expressed_transcripts_bed_parameters,
      params.lens$neos$indels_to_neos$samtools_faidx_fetch_somatic_folder_parameters)

    // Self-antigens to peptides
    selfs_to_neos(
      alns_to_transcript_counts.out.quants,
      germ_and_som_vars_to_tumor_vars.out.tumor_vars,
      params.lens$neos$selfs_to_neos$gtf,
      params.lens$neos$selfs_to_neos$dna_ref,
      params.lens$neos$selfs_to_neos$cta_self_gene_list,
      params.lens$neos$selfs_to_neos$samtools_index_parameters,
      params.lens$neos$selfs_to_neos$lenstools_filter_expressed_self_parameters,
      params.lens$neos$selfs_to_neos$lenstools_get_expressed_self_bed_parameters,
      params.lens$neos$selfs_to_neos$samtools_faidx_fetch_parameters,
      params.lens$neos$selfs_to_neos$bcftools_index_parameters)

    // ERVs to peptides
    ervs_to_neos(
      params.lens$neos$ervs_to_neos$dna_ref,
      samtools_index_rna.out.bams_and_bais,
      samtools_index_txome.out.bams_and_bais,
      alns_to_transcript_counts.out.quants,
      params.lens$neos$ervs_to_neos$geve_general_ref,
      params.lens$neos$ervs_to_neos$lenstools_get_expressed_ervs_bed_parameters,
      params.lens$neos$ervs_to_neos$lenstools_make_erv_peptides_parameters,
      params.lens$neos$ervs_to_neos$lenstools_filter_expressed_ervs_parameters,
      params.lens$neos$ervs_to_neos$lenstools_filter_ervs_by_rna_coverage_parameters,
      params.lens$neos$ervs_to_neos$normal_control_quant,
      params.lens$neos$ervs_to_neos$tpm_threshold)

    // Viruses to peptides
    viruses_to_neos(
      unaligned_fqs_to_virdetect_cds_counts.out.viral_cds_counts,
      unaligned_fqs_to_virdetect_cds_counts.out.viral_cds_alns,
      params.lens$neos$viruses_to_neos$viral_cds_ref,
      params.lens$neos$viruses_to_neos$lenstools_filter_expressed_viruses_parameters,
      params.lens$neos$viruses_to_neos$lenstools_filter_viruses_by_rna_coverage_parameters,
      params.lens$neos$viruses_to_neos$lenstools_get_expressed_viral_bed_parameters,
      params.lens$neos$viruses_to_neos$lenstools_make_viral_peptides_parameters)

    // Splice variants to peptides
    // Currently handled by alns_to_splice_variants due to NeoSplice being an
    // all inclusive tool. This step will be needed in the future when using
    // other splice variant tools.

    // Fusion variants to peptides
    fusions_to_neos(
      procd_fqs_to_fusions.out.coding_effect_fusions,
      make_phased_germline_vars.out.phased_vcfs,
      params.lens$neos$fusions_to_neos$gtf,
      params.lens$neos$fusions_to_neos$dna_ref,
      params.lens$neos$fusions_to_neos$bedtools_index_phased_germline_parameters,
      params.lens$neos$fusions_to_neos$lenstools_get_fusion_transcripts_bed_parameters,
      params.lens$neos$fusions_to_neos$samtools_faidx_fetch_parameters)

    // Combining peptides from all antigen sources
    selfs_to_neos.out.self_antigen_c1_peptides.map{ [it[0], it[2], it[3]] }
      .join(ervs_to_neos.out.erv_c1_peptides.map{ [it[0], it[2], it[3]] }, by: [0,1], remainder: true)
      .join(snvs_to_neos.out.som_var_c1_peptides.map{ [it[0], it[3], it[4]] }, by: [0,1], remainder: true)
      .join(indels_to_neos.out.som_var_c1_peptides.map{ [it[0], it[3], it[4]] }, by: [0,1], remainder: true)
      .join(fusions_to_neos.out.fusion_c1_peptides.map{ [it[0], it[3], it[4]] }, by: [0,1], remainder: true)
      .join(viruses_to_neos.out.viral_c1_peptides.map{ [it[0], it[2], it[3]] }, by: [0,1], remainder: true)
      .join(alns_to_splice_variants.out.splice_c1_peptides.map{ [it[0], it[3], it[4]] }, by: [0,1], remainder: true)
      .map{ [*it.asList().minus(null)] }
      .map{ [it[0], it[1], it[2..-1]] }
      .set{ init_joint_peptides }
    combine_peptide_fastas(
      init_joint_peptides)
      .map{ [it[0], 'NA', it[1], it[2]] }
      .set{ combined_peptides }

    // Combining NT sequence from ERVs, Viruses, and CTA/Self-antigens for quantification
    selfs_to_neos.out.self_antigen_c1_nts.map{ [it[0], it[2], it[3]] }
      .join(ervs_to_neos.out.erv_c1_nts.map{ [it[0], it[2], it[3]] }, by: [0,1], remainder: true)
      .join(viruses_to_neos.out.viral_c1_nts.map{ [it[0], it[2], it[3]] }, by: [0,1], remainder: true)
      .map{ [*it.asList().minus(null)] }
      .map{ [it[0], it[1], it[2..-1]] }
      .set{ init_joint_nts }
    combine_nt_fastas(
      init_joint_nts)
      .map{ [it[0], 'NA', it[1], it[2]] }
      .set{ combined_nts }


    // Summarize pMHCs
    peps_and_alleles_to_antigen_stats(
      combined_peptides,
      all_pat_mhc_alleles,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$antigen_tool,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$antigen_tool_parameters,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$antigen_tool_ref_dir,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$species,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$peptide_lengths)

    aggregate_pmhc_summaries(
      peps_and_alleles_to_antigen_stats.out.all_antigen_outputs)

    aggregate_pmhc_summaries.out.pmhc_aggs_tsv
      .join(combine_peptide_fastas.out.peptide_fastas, by: [0,1])
      .set{ pmhc_summs_and_fastas }

    lenstools_annotate_pmhcs(
      pmhc_summs_and_fastas,
      params.lens$lenstools$lenstools_annotate_pmhcs$binding_affinity_threshold)



    erv_read_counts = Channel.empty()
    cta_read_counts = Channel.empty()
    snv_read_counts = Channel.empty()
    indel_read_counts = Channel.empty()
    fusion_read_counts = Channel.empty()
    splice_read_counts = Channel.empty()
    viral_read_counts = Channel.empty()
*/

/*
    // Get CTA/self-antigen and ERV RNA read support
    lenstools_get_erv_and_cta_peptide_read_count(
      samtools_index_txome.out.bams_and_bais,
      combined_nts,
      lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs)
    lenstools_get_erv_and_cta_peptide_read_count.out.erv_peptide_read_counts
      .set{ erv_read_counts }
    lenstools_get_erv_and_cta_peptide_read_count.out.self_peptide_read_counts
      .set{ cta_read_counts }
*/

/*
    // Get viral RNA read support
    viruses_to_neos.out.viral_bams_and_bais.map{ [it[0], it[2], it[3], it[4]] }
      .join(viruses_to_neos.out.viral_consensus_fastas.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs, by: [0, 1])
      .set{ viral_peptide_read_count_inputs }
    lenstools_get_viral_peptide_read_count(
      viral_peptide_read_count_inputs)
    lenstools_get_viral_peptide_read_count.out.peptide_read_counts
      .set{ viral_read_counts }
*/

/*
    // Get SNV RNA read support
    snvs_to_neos.out.som_var_c1_nts.map{ [it[0], it[3], it[4]] }
      .join(lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs, by: [0, 1])
      .join(samtools_index_txome.out.bams_and_bais, by: 0)
      .map{ [it[0], it[1], it[6], it[7], it[2], it[3]] }
      .set{ snv_peptide_read_count_inputs }
    lenstools_get_snv_peptide_read_count(
      snv_peptide_read_count_inputs,
      params.lens$lenstools$lenstools_get_snv_peptide_read_count$gtf)
    lenstools_get_snv_peptide_read_count.out.peptide_read_counts
      .set{ snv_read_counts }

    // Get InDel RNA read support
    indels_to_neos.out.som_var_c1_nts.map{ [it[0], it[3], it[4]] }
      .join(lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs, by: [0, 1])
      .join(samtools_index_rna.out.bams_and_bais, by: 0)
      .map{ [it[0], it[1], it[6], it[7], it[2], it[3]] }
      .set{ indel_peptide_read_count_inputs }
    lenstools_get_indel_peptide_read_count(
      indel_peptide_read_count_inputs,
      params.lens$lenstools$lenstools_get_indel_peptide_read_count$gtf)
    lenstools_get_indel_peptide_read_count.out.peptide_read_counts
      .set{ indel_read_counts }

    // Get splice RNA read support
    lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs
      .join(samtools_index_rna.out.bams_and_bais.map{ [it[0], it[2], it[3], it[4]] }, by: [0, 1])
      .set{ splice_peptide_read_count_inputs }
    lenstools_get_splice_peptide_read_count(
      splice_peptide_read_count_inputs)
    lenstools_get_splice_peptide_read_count.out.peptide_read_counts
      .set{ splice_read_counts }

    // Get fusion RNA read support
    fusions_to_extracted_reads(
      procd_fqs_to_fusions.out.fusions)
    manifest_to_rna_alns.out.procd_fqs.filter{ it[1] =~ 'ar-' }
      .join(fusions_to_extracted_reads.out.fusion_read_names, by: [0, 1, 2])
      .set{ tumor_reads_and_fusion_read_names }
    seqtk_subseq(
      tumor_reads_and_fusion_read_names,
      '.fusion_reads',
      params.lens$seqtk_subseq_parameters)
    procd_fqs_to_fusions.out.fusions
      .join(seqtk_subseq.out.extracted_fqs, by: [0, 1, 2])
      .join(fusions_to_neos.out.fusion_c1_nts, by: [0, 1, 2])
      .map{ [it[0], it[2], it[3], it[4], it[5]] }
      .join(lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs, by: [0, 1])
      .set{ fusion_peptide_read_count_inputs }
    lenstools_get_fusion_peptide_read_count(
      fusion_peptide_read_count_inputs)
    lenstools_get_fusion_peptide_read_count.out.peptide_read_counts
      .set{ fusion_read_counts }


    Channel.empty()
     .concat(snv_read_counts)
     .concat(indel_read_counts)
     .concat(splice_read_counts)
     .concat(fusion_read_counts)
     .concat(erv_read_counts)
     .concat(cta_read_counts)
     .concat(viral_read_counts)
     .groupTuple(by: [0, 1])
     .set{ pan_source_counts_by_pat }

   lenstools_combine_read_counts(
     pan_source_counts_by_pat)


    // Need to capture read counts here.
    // Agretopicity
    calculate_agretopicity(
      lenstools_combine_read_counts.out.pmhcs_with_read_counts,
      params.lens$immuno$calculate_agretopicity$blastp_db_dir,
      params.lens$immuno$peps_and_alleles_to_antigen_stats$species,
      all_pat_mhc_alleles,
      params.lens$immuno$calculate_agretopicity$peps_and_alleles_to_antigen_stats$antigen_tool,
      params.lens$immuno$calculate_agretopicity$peps_and_alleles_to_antigen_stats$antigen_tool_parameters,
      params.lens$immuno$calculate_agretopicity$peps_and_alleles_to_antigen_stats$antigen_tool_ref_dir,
      params.lens$immuno$calculate_agretopicity$peps_and_alleles_to_antigen_stats$peptide_lengths)

    calculate_agretopicity.out.pmhcs_with_agretos
      .join(cnas_and_vcfs_to_ccfs.out.ccfs.map{ [it[0], it[3], it[4]] }, by: [0, 1])
      .set{ annot_ccfs_input }
    annotate_ccfs(
      annot_ccfs_input,
      '.ccfs',
      '-a variant_pos -b mutation_id -c cluster_id,cellular_prevalence,cluster_assignment_prob')


    // Make antigen-specific IGV files
    lenstools_make_lens_bed(
        lenstools_annotate_pmhcs.out.high_aff_annotated_pmhcs)

    samtools_index_dna(
      alns_to_dna_procd_alns.out.procd_bams,
      '')

    samtools_index_dna.out.bams_and_bais.filter{ it[1] =~ 'nd-' }.set{ norm_dna_bams_bais }
    samtools_index_dna.out.bams_and_bais.filter{ it[1] =~ 'ad-' }.set{ tumor_dna_bams_bais }

    norm_dna_bams_bais
      .join(tumor_dna_bams_bais, by: [0, 2])
      .set{ norm_tumor_dna_bams_and_bais }

    norm_tumor_dna_bams_and_bais
      .join(alns_to_rna_procd_alns.out.procd_bams_and_bais
        .filter{ it[1] =~ 'ar-' }
        .map{ [it[0], it[2], it[1], it[3], it[4]] }, by: [0, 1])
      .set{ all_pat_bams }

    all_pat_bams
      .join(lenstools_make_lens_bed.out.lens_bed, by: [0, 1])
      .set{ all_pat_bams_with_beds }

    bam_subsetter(
      all_pat_bams_with_beds)

    igv_snapshot_automator(
      bam_subsetter.out.subsetted_bams_w_bed)

    lenstools_make_lens_report(
      annotate_ccfs.out.annoted_files,
      params.lens_out_dir)
*/
    }
