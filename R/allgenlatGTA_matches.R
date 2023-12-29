#' Collects gen/lat/GTA classification pattern-match information
#'
#' Collects pattern information associated with all contigs predicted as containing gen/lat/GTA transduction events based on the results from the pattern_matcher function.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @keywords internal
allgenlatGTA_matches <- function(best_match_list){
  A<-1
  genlatGTAclassification_list <- list()
  for (i in seq(1,length(best_match_list),1)){
    classification <-  best_match_list[[i]][[8]]
    if(classification !="Gen/Lat/GTA") next
    genlatGTAclassification_list[[A]] <- best_match_list[[i]]
    A <- A+1
  }
  return(genlatGTAclassification_list)
}
