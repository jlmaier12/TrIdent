#' Summarizes pattern-matches into a table
#'
#' Summarizes the classications made in the pattern_matcher function based on which pattern match achieved the lowest mean absolute difference for each contig. Outputs results in a table.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @keywords internal
contig_classification_summary <- function(best_match_list){
  ref_name <- rep(NA, length(best_match_list))
  classifications <- rep(NA, length(best_match_list))
  NormMatchScore <- rep(NA, length(best_match_list))
  lapply(1:length(best_match_list), function(i) {
    NormMatchScore[i] <<- best_match_list[[i]][[10]]
    ref_name[i] <<- best_match_list[[i]][[9]]
    classifications[i] <<- best_match_list[[i]][[8]]
  })
  classification_summary <- cbind.data.frame(ref_name, classifications, NormMatchScore)
  return(classification_summary)
}
