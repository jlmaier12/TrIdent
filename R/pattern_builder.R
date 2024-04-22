#' Pattern-builder
#'
#' Builds the pattern (vector) associated with the 'best pattern-match' for each contig classified as prophage-like, gen/lat/GTA, or none with high VLP:WC read coverage
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param transductionclassifications_wlowratios A list containing pattern match information associated with all contigs classified as prophage-like, gen/lat/GTA, or none with high VLP-fraction:whole-community read coverage ratios. Generated with the alltransduction_withlowratios_summarylist
#' @param i The index position associated with the current contig's best pattern-match information
#' @param classification The classification of the contig assigned by the main TrIdent_classifier function
#' @keywords internal
pattern_builder <- function(viral_subset, transductionclassifications_wlowratios, classification, i){
  min_read_cov <- transductionclassifications_wlowratios[[i]][[2]]
  max_read_cov <- transductionclassifications_wlowratios[[i]][[3]]
  gen_covsteps <- transductionclassifications_wlowratios[[i]][[4]]
  start_pos <- transductionclassifications_wlowratios[[i]][[5]]
  end_pos <- transductionclassifications_wlowratios[[i]][[6]]
  if (classification=="Prophage-like"){
    if (start_pos==1) {
      pattern <- c(rep(max_read_cov,end_pos), rep(min_read_cov, (nrow(viral_subset)-end_pos)))
    } else if (end_pos == nrow(viral_subset)){
      pattern <- c(rep(min_read_cov, start_pos), rep(max_read_cov, (nrow(viral_subset)-start_pos)))
    } else{
      match_region <- end_pos-start_pos
      pattern <- c(rep(min_read_cov, start_pos), rep(max_read_cov, match_region), rep(min_read_cov, (nrow(viral_subset)-(match_region+start_pos))))
    }
  } else if (classification=="Sloping") {
    if (start_pos==1 & end_pos == nrow(viral_subset) & gen_covsteps < 0){
      pattern <- seq(max_read_cov, min_read_cov, gen_covsteps)
    } else if (start_pos==1 & end_pos == nrow(viral_subset) & gen_covsteps > 0){
      pattern <- seq(min_read_cov, max_read_cov, gen_covsteps)
    } else if (start_pos != 1){
      pattern <- c(rep(min(viral_subset[,2]), start_pos-1), seq(max_read_cov, min_read_cov, gen_covsteps))
    } else if (end_pos != nrow(viral_subset)) {
      pattern <- c(seq(min_read_cov, max_read_cov, gen_covsteps),rep(min(viral_subset[,2]), (nrow(viral_subset)-end_pos)))
    }
  } else if (classification=="HighCoverageNoPattern") {
    pattern <- rep(min_read_cov, max_read_cov) #for nones, min read cov is actually mean read cov and max read cov is the contig length
  }
  return(pattern)
}
