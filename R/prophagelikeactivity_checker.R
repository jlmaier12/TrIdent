#' Prophage-like activity/abundance checker
#'
#' Determines whether a detected prophage/prophage-like element is active/highly abundant based on elevated (>1.2x the average read coverage of the non-prophage-like region) read coverage of the prophage-like region in the whole-community reads.
#'
#' @param classificationsummary Classification summary table
#' @param prophageclassifications A list containing pattern match information associated with all contigs classified as prophage-like. Generated with the allprophagelike_matches
#' @param phageread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
prophagelikeactivity_checker <- function(classificationsummary, prophageclassifications, phageread_dataset, microbialread_dataset, windowsize){
  classificationsummary$active_prophage <- rep(NA, nrow(classificationsummary))
  classificationsummary$elevation_ratio <- rep(NA, nrow(classificationsummary))
  if(length(prophageclassifications)==0) {
    return(classificationsummary)
  }
  lapply(1:length(prophageclassifications), function(i) {
    viral_subset <- phageread_dataset[which(phageread_dataset[,1] == prophageclassifications[[i]][[9]]),]
    viral_subset <- windowsize_func(viral_subset,windowsize)
    start_pos <- prophageclassifications[[i]][[5]]
    end_pos <- prophageclassifications[[i]][[6]]
    ref_name <- prophageclassifications[[i]][[9]]
    match_length <- abs(end_pos-start_pos) * windowsize
    nonmatch_lengthbp <- (nrow(viral_subset)*windowsize)-match_length
    if (nonmatch_lengthbp < 20000) {
      classificationsummary[which(classificationsummary[,1]==ref_name),6] <<- "CBD" #cant be determined
      classificationsummary[which(classificationsummary[,1]==ref_name),7] <<- "CBD"
    } else {
      microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == prophageclassifications[[i]][[9]]),]
      microbial_subset <- windowsize_func(microbial_subset,windowsize)
      prophage_region <- microbial_subset[c(start_pos:end_pos),2]
      nonprophage_region <- microbial_subset[c(c(0:start_pos),c(end_pos:nrow(microbial_subset))),2]
      ratio <- round(mean(prophage_region)/mean(nonprophage_region), digits=4)
      classificationsummary[which(classificationsummary[,1]==ref_name),7] <<- ratio
      if(ratio > 1.3) {
        classificationsummary[which(classificationsummary[,1]==ref_name),6] <<- "YES"
      } else if (ratio < 0.5){
        classificationsummary[which(classificationsummary[,1]==ref_name),6] <<- "MIXED"
      } else {
        classificationsummary[which(classificationsummary[,1]==ref_name),6] <<- "NO"
      }
    }
  })
  return(classificationsummary)
}
