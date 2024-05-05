#' @export
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
    print("max colsums is 3")
    # get cols with that colsum
    colsum_3 <- which(colSums(comparison) == 3)
    tool_name <- colnames(colsum_3)[1]
    output <- calls[[tool_name]]
    # there is a tie, max colsum is 2
  } else {
    print("2 callers give same result")
    # check if total sum is 8 - if so, we have two equal ties - go for optitype
    if (sum(comparison) == 8) {
      output <- calls[["optitype"]]
    } else {
      # otherwise, get the name of one of the tools in the tie and go for the call
      colsum_2 <- which(colSums(comparison) == 2)
      tool_name <- colnames(comparison[,colsum_2])[1]
      print("tool name to get call for")
      print(tool_name)
      output <- calls[[tool_name]]
    }
  }
  return(output)
}

#' @export
majority_vote2 <- function(comparison, calls) {
  stopifnot(nrow(comparison) == 4)
  stopifnot(length(calls) == 4)
  stopifnot(all(c("optitype", "polysolver", "kourami", "hlala") %in% rownames(comparison)))
  # check if there are any NAs
  print(calls)
  if (!anyNA(calls, recursive = TRUE)){
	  print("no NAs")
	  # all different, go with optitype (best tool for all MHC Class I alleles)
	  if (max(colSums(comparison)) == 1) {
	    output <- calls[["optitype"]]
	    # all agree, doesn't matter which one we go with
	  } else if (max(colSums(comparison)) == 4) {
	    output <- calls[["optitype"]]
	  } else if (max(colSums(comparison)) == 3) {
	    print("max colsums is 3")
	    # get cols with that colsum
	    colsum_3 <- which(colSums(comparison) == 3)
	    tool_name <- colnames(colsum_3)[1]
	    output <- calls[[tool_name]]
	    # there is a tie, max colsum is 2
	  } else {
	    print("2 callers give same result")
	    # check if total sum is 8 - if so, we have two equal ties - go for optitype
	    if (sum(comparison) == 8) {
	      output <- calls[["optitype"]]
	    } else {
	      # otherwise, get the name of one of the tools in the tie and go for the call
	      colsum_2 <- which(colSums(comparison) == 2)
	      tool_name <- colnames(comparison[,colsum_2])[1]
	      print("tool name to get call for")
	      print(tool_name)
	      output <- calls[[tool_name]]
	    }
	 }
  } else {
    print("at least one element NA")
    output <- c(NA, NA)
  }
  return(output)
}

