#' No transduction pattern match
#'
#' Assess whether a contig does not have a read coverage pattern associated with a transduction event. A horizontal line at the mean coverage should be an optimal match if the contig read coverage displays no transduction patterns
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @keywords internal
notransduction_pattern <- function (viral_subset) {
  pattern <- rep(mean(viral_subset[,2]),nrow(viral_subset))
  Cov_values_contig <- viral_subset[,2]
  diff <- mean(abs(Cov_values_contig - pattern))
  best_match_info <- list(diff, mean(viral_subset[,2]), nrow(viral_subset), "NA", "NA", "NA", "NA","InsufficientCoverage")
  return(best_match_info)
}
