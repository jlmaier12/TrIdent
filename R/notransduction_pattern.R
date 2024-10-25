#' No transduction pattern match
#'
#' Assess whether a contig does not have a read coverage pattern associated with a transduction event. A horizontal line at the mean coverage should be an optimal match if the contig read coverage displays no transduction patterns
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @keywords internal
notransduction_pattern <- function (viral_subset) {
  pattern1 <- rep(median(viral_subset[,2]),nrow(viral_subset))
  pattern2 <- rep(mean(viral_subset[,2]),nrow(viral_subset))
  Cov_values_contig <- viral_subset[,2]
  diff1 <- mean(abs(Cov_values_contig - pattern1))
  diff2 <- mean(abs(Cov_values_contig - pattern2))
  diff <- ifelse(diff1 < diff2, diff1, diff2)
  best_match_info <- list(diff, mean(viral_subset[,2]), nrow(viral_subset), "NA", "NA", "NA", "NA","InsufficientCoverage")
  return(best_match_info)
}
