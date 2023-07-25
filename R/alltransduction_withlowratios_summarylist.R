#' Collects Pattern-match information associated with prophage-like and gen/lat/GTA classifications and contigs with low whole-community:VLP-fraction read coverage ratios
#'
#' Collects Pattern information associated with all contigs classified as prophage-like or gen/lat/GTA based on the results from the pattern_matcher function and all contigs with 'none' classifications that have low whole-community:VLP-fraction read coverage ratios.
#'
#' @param best_match_list Classifications made with pattern_matcher function. Classifications are stored as the first item in the best_match_list.
#' @param summarytable_wratios Classification summary table with whole-community:VLP-fraction read coverage ratios calculated.
#' @keywords internal
alltransduction_withlowratios_summarylist <- function(best_match_list, summarytable_wratios){
  A<-1
  transductionclassification_withlowratios_list <- list()
  for (i in seq(1,length(best_match_list),1)){
    ref_name <- best_match_list[[i]][[8]]
    classification <-  summarytable_wratios[which(summarytable_wratios[,1]==ref_name),2]
    if(classification=="None") next
    transductionclassification_withlowratios_list[[A]] <- best_match_list[[i]]
    A <- A+1
  }
  return(transductionclassification_withlowratios_list)
}
