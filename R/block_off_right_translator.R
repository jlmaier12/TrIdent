#' Block off-right pattern translator
#'
#' Build and translate a block pattern going off the right side of the contig. Stop translating when the pattern left on the contig is 10,000bp. Translate the pattern 2000bp at a time.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
block_off_right_translator <- function (viral_subset, windowsize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  quarter_read_cov <- abs(max_read_cov-min_read_cov)/4
  bottomtotop_read_cov <- abs(max_read_cov - quarter_read_cov)/10
  Cov_values_contig <- viral_subset[,2]
  if(min_read_cov>(max_read_cov/10)) {
    min_read_cov <- 0
  }
  startingcoverages <- seq(quarter_read_cov,  max_read_cov, bottomtotop_read_cov)
  pattern <- c(rep(min_read_cov, 10000/windowsize), rep(startingcoverages[1], nrow(viral_subset)-(10000/windowsize)))
  diff <- mean(abs(Cov_values_contig - pattern))
  start_pos <- (which(pattern == max(pattern))[1])
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], "NA", start_pos, length(pattern))
  for(cov in startingcoverages) {
    pattern <- c(rep(min_read_cov, 10000/windowsize), rep(cov, nrow(viral_subset)-(10000/windowsize)))
    repeat {
      if (length(which(pattern==cov)) < 10000/windowsize) break
      if (diff < best_match_info[[1]]){
        best_match_info <- list(diff, min_read_cov, cov, "NA", start_pos, length(pattern))
      }
      pattern <- c(rep(min_read_cov,(2000/windowsize)),pattern[-c(((length(pattern))-((2000/windowsize)-1)):length(pattern))]) #variable, removing 2000bp at a time
      diff <- mean(abs(Cov_values_contig - pattern))
      start_pos <- (which(pattern == max(pattern))[1])
    }
  }
  new_pattern_coverages <- blockheight_optim(best_match_info, startingcoverages, bottomtotop_read_cov)
  for(newcov in new_pattern_coverages) {
    pattern <- c(rep(min_read_cov, 10000/windowsize), rep(newcov, nrow(viral_subset)-(10000/windowsize)))
    repeat {
      if (length(which(pattern==newcov)) < 10000/windowsize) break
      if (diff < best_match_info[[1]]){
        best_match_info <- list(diff, min_read_cov, cov, "NA", start_pos, length(pattern))
      }
      pattern <- c(rep(min_read_cov,(2000/windowsize)),pattern[-c(((length(pattern))-((2000/windowsize)-1)):length(pattern))]) #variable, removing 2000bp at a time
      diff <- mean(abs(Cov_values_contig - pattern))
      start_pos <- (which(pattern == max(pattern))[1])
    }
  }
  best_match_results <- append(best_match_info, "Prophage-like")
  return(best_match_results)
}
