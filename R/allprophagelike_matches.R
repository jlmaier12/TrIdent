#' Collects prophage-like classification pattern-match information
#'
#' Collects pattern information associated with all contigs predicted as containing prophages based on the results from the pattern_matcher function.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @keywords internal
allprophagelike_matches <- function(best_match_list){
  A<-1
  prophagelikeclassification_list <- list()
  for (i in seq(1,length(best_match_list),1)){
    classification <-  best_match_list[[i]][[8]]
    if(classification !="Prophage-like") next
    prophagelikeclassification_list[[A]] <- best_match_list[[i]]
    A <- A+1
  }
  return(prophagelikeclassification_list)
}
