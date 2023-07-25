#' Pattern-match size calculator
#'
#' Calculate the size, in base pairs, of the matching region for prophage-like and gen/lat/GTA patterns
#'
#' @param classificationsummary Classification summary table
#' @param transductionclassification_list A list containing pattern match information associated with all contigs classified as prophage-like or gen/lat/GTA. Generated with the alltransduction_events_summarylist
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
matchsize_checker <- function(classificationsummary, transductionclassification_list, windowsize){
  classificationsummary <- as.data.frame(classificationsummary)
  classificationsummary$match_size <- rep(NA, nrow(classificationsummary))
  for (i in seq(1,length(transductionclassification_list),1)) {
    ref_name <- transductionclassification_list[[i]][[8]]
    start_pos <- transductionclassification_list[[i]][[5]]
    end_pos <- transductionclassification_list[[i]][[6]]
    classificationsummary[which(classificationsummary[,1]==ref_name),3] <- (end_pos-start_pos) *windowsize
  }
  return(classificationsummary)
}
