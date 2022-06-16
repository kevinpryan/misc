# some useful functions for very specific use cases

# function to parse EGAX file, extracting Erx ID and library construction protocol
extract_erx_id_lib_construction <- function(path_to_xml) {
  require(xml2)
  require(dplyr)
  file_in <- read_xml(path_to_xml) %>% as_list()
  lib_prep_kit <- file_in$EXPERIMENT_SET$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR$LIBRARY_CONSTRUCTION_PROTOCOL[[1]]
  erx_no <- file_in$EXPERIMENT_SET$EXPERIMENT$IDENTIFIERS$PRIMARY_ID[[1]]
  output <- c(erx_no,lib_prep_kit)
  return(output)
}

# function to parse runs/EGAR file, extracting ERX ID and library construction protocol
extract_erx_id_sample_no <- function(path_to_xml) {
  require(xml2)
  require(dplyr)
  require(stringr)
  file_in <- read_xml(path_to_xml) %>% as_list()
  erx_no <- file_in$RUN_SET$RUN$EXPERIMENT_REF$IDENTIFIERS$PRIMARY_ID
  file_attributes_r1 <- attributes(file_in$RUN_SET$RUN$DATA_BLOCK$FILES$FILE)
  name_file_r1 <- file_attributes_r1$filename
  samplename_r1 <- str_split_fixed(name_file_r1, pattern = "/", n = 2)[,2]
  samplename_r2 <- sub(pattern = ".R1.", replacement = ".R2.", x = samplename_r1)
  output <- c(erx_no, samplename_r1, samplename_r2)
  return(output)
}