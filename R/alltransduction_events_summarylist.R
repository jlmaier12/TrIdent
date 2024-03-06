#' Collects prophage-like and gen/lat/GTA transduction classification pattern-match information
#'
#' Collects pattern information associated with all contigs classified as prophage-like or gen/lat/GTA based on the results from the pattern_matcher function.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @keywords internal
alltransduction_events_summarylist <- function(best_match_list){
  A<-1
  transductionclassification_list <- list()
  lapply(1:length(best_match_list), function(i) {
    classification <-  best_match_list[[i]][[8]]
    if(classification=="None") return(NULL)
    transductionclassification_list[[A]] <<- best_match_list[[i]]
    A <<- A+1
  })
  return(transductionclassification_list)
}
