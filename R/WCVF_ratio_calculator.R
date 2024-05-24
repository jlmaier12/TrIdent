#' Whole-community:VLP-fraction read coverage ratio calculator
#'
#' Calculate the whole-community:VLP-fraction read coverage ratio using the median coverage for all contigs. If the ratio is less than 2 (i.e VLP-fraction read coverage is on average at least half the whole-community read coverage), then mark the contig as having a low whole-community (WC) to VLP-fraction (VLPF) read coverage ratio.
#'
#' @param classificationsummary Classification summary table
#' @param microbialread_dataset A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param phageread_dataset  A table containing contig names, coverages averaged over 100bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @keywords internal
WCVF_ratio_calculator <- function(classificationsummary, microbialread_dataset, phageread_dataset){
  None_indexes <- which(classificationsummary[,2]=="InsufficientCoverage")
  lapply(1:length(None_indexes), function(p) {
    i<-None_indexes[[p]]
    ref_name <- classificationsummary[i,1]
    viral_subset <- phageread_dataset[which(phageread_dataset[,1] == ref_name),]
    viral_subset[is.nan.data.frame(viral_subset)] <- 0
    microbial_subset <- microbialread_dataset[which(microbialread_dataset[,1] == ref_name),]
    microbial_subset[is.nan.data.frame(microbial_subset)] <- 0
    WC_to_VLPF_ratio <- median(microbial_subset[,2])/ median(viral_subset[,2])
    if(WC_to_VLPF_ratio < 2) {
      classificationsummary[i,2] <<- "HighCoverageNoPattern"
    } else {
      classificationsummary[i,2] <<- "InsufficientCoverage"
    }
  })
  return(classificationsummary)
}
