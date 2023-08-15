#!/usr/bin/env nextflow
params.lens$somatic$combine_strategy = "intersect"
params.lens$somatic$alns_to_som_vars$som_var_caller = "strelka2,mutect2,abra2"
//params.lens$somatic$alns_to_som_vars$som_var_caller = "strelka2"
    // Combine variants called from different variant callers using either intersection or union
    // if neither strategy is chosen, just set joint_vars to the output of htslib_bgzip_somatic
    if (params.lens$somatic$combine_strategy =~ /intersect/) {
      // Somatic variant intersectioning
      if (params.lens$somatic$alns_to_som_vars$som_var_caller =~ /,/){
        println "Intersecting variants from variant callers"
        /*
        som_vars_to_isecd_som_vars(
          htslib_bgzip_somatic.out.bgzip_files,
          params.lens$somatic$som_vars_to_isecd_som_vars$vcf_isecing_tool,
          params.lens$somatic$som_vars_to_isecd_som_vars$vcf_isecing_tool_parameters)
        som_vars_to_isecd_som_vars.out.isecd_som_vars
          .set{ joint_vars }
        */
     } else {
        error "ERROR: params.lens somatic combine_strategy is intersect, but multiple variant callers not used."
      }
     }
     else if (params.lens$somatic$combine_strategy =~ /merge|union/) {
       if (params.lens$somatic$alns_to_som_vars$som_var_caller =~ /,/){
        println "merging/getting union of variants from variant callers"
        /*
        som_vars_to_union_som_vars(
          htslib_bgzip_somatic.out.bgzip_files,
          params.lens$somatic$som_vars_to_union_som_vars$vcf_merging_tool,
          params.lens$somatic$som_vars_to_union_som_vars$vcf_merging_tool_parameters)
        som_vars_to_union_som_vars.out.union_som_vars
          .set{ joint_vars }
        */
        } else {
          error "ERROR: params.lens somatic combine_strategy is merge/union, but multiple variant callers not included in params.lens somatic alns_to_som_vars som_var_caller"
        }
    } else {
        if (params.lens$somatic$alns_to_som_vars$som_var_caller =~ /,/) {
           error "ERROR: 'params.lens somatic combine_strategy' chosen is either not valid or is empty, and multiple params.lens somatic alns_to_som_vars som_var_callers have been chosen. Please use one of intersect, merge or union or choose one variant caller only."
       } else {
           println "No combination strategy specified, continuing with the one variant caller used"
           //som_vars_to_normd_som_vars.out.normd_som_vars.set{ joint_vars }
       }
    }


