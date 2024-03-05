# function that takes in a dataframe in which the rownames are your versioned ensemble gene ids (e.g. ENSG12345.12) and returns a dataframe in which the rownames are hgnc symbols (e.g. COL6A3)
# if the new name is not found (NA) or is duplicated (i.e. multiple HGNC symbols for the same ENSG symbol), the ENSG symbol is retained in the output
convert_ensg_version_to_hgnc_df_rownames <- function(df, ensdb = EnsDb.Hsapiens.v86){
  library(AnnotationDbi)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  genes_base <- str_split_fixed(string = rownames(df), pattern = "\\.", n = 2)[,1]
  newnames_original <- suppressWarnings(mapIds(EnsDb.Hsapiens.v86,
    keys = genes_base,
    column = 'SYMBOL',
    keytype = 'GENEID'))

# keep ensg version of newnames is na or is duplicated
newnames <- ifelse(is.na(newnames_original) | duplicated(newnames_original),
    rownames(df), newnames_original)
rownames(df) <- newnames
return(df)
}
