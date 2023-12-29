#' Summarize slopes for Gen/Lat/GTA classifications
#'
#' Add slope information for Gen/Lat/GTA classifications to summary table
#'
#' @param classificationsummary Classification summary table
#' @param genlatGTAclassifications A list containing pattern match information associated with all contigs classified as gen/lat/GTA. Generated with the allgenlatGTA_matches
#' @keywords internal
slope_summary <- function(classificationsummary, genlatGTAclassifications){
  classificationsummary$slope <- rep(NA, nrow(classificationsummary))
  if(length(genlatGTAclassifications)==0) {
    return(classificationsummary)
  }
  for (i in seq(1,length(genlatGTAclassifications),1)) {
    slope <- round(genlatGTAclassifications[[i]][[7]], digits=4)
    ref_name <- genlatGTAclassifications[[i]][[9]]
    classificationsummary[which(classificationsummary[,1]==ref_name),8] <- slope
  }
  return(classificationsummary)
}
