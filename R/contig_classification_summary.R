#' Summarizes pattern-matches into a table
#'
#' Summarizes the classications made in the pattern_matcher function based on which pattern match achieved the lowest mean absolute difference for each contig. Outputs results in a table.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @keywords internal
contig_classification_summary <- function(best_match_list){
  ref_name <- rep(NA, length(best_match_list))
  classifications <- rep(NA, length(best_match_list))
  for (i in seq(1,length(best_match_list),1)){
    ref_name[i] <- best_match_list[[i]][[8]]
    classifications[i] <- best_match_list[[i]][[7]]
  }
  classification_summary <- cbind(ref_name, classifications)
  return(classification_summary)
}
