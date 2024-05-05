majority_vote <- function(comparison, calls) {
  stopifnot(nrow(comparison) == 4)
  stopifnot(length(calls) == 4)
  stopifnot(all(c("optitype", "polysolver", "kourami", "hlala") %in% rownames(comparison)))
  # all different, go with optitype (best tool for all MHC Class I alleles)
  if (max(colSums(comparison)) == 1) {
    output <- calls[["optitype"]]
    # all agree, doesn't matter which one we go with
  } else if (max(colSums(comparison)) == 4) {
    output <- calls[["optitype"]]
  } else if (max(colSums(comparison)) == 3) {
    # get cols with that colsum
    colsum_3 <- which(colSums(comparison) == 3)
    tool_name <- colnames(colsum_3)[1]
    output <- calls[[tool_name]]
    # there is a tie, max colsum is 2
  } else {
    # check if total sum is 8 - if so, we have two equal ties - go for optitype
    if (sum(comparison) == 8) {
      output <- calls[["optitype"]]
    } else {
      # otherwise, get the name of one of the tools in the tie and go for the call
      colsum_2 <- which(colSums(comparison) == 2)
      tool_name <- colnames(colsum_2)[1]
      output <- calls[[tool_name]]
    }
  }
  return(output)
}
