#' Builds pattern for full block-like pattern
#'
#' Build the full block pattern and provide it as input to the blockpattern_translator to be translated across a contig. The pattern is made smaller length-wise 2000bp at a time, and the pattern stops decreasing in length once it reaches 10000bp.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param blocksize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @keywords internal
full_blockpattern_builder <- function (viral_subset, windowsize, blocksize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  quarter_read_cov <- abs(max_read_cov-min_read_cov)/4
  bottomtotop_read_cov <- abs(max_read_cov - (min_read_cov+quarter_read_cov))/10
  Cov_values_contig <- viral_subset[,2]
  if(min_read_cov>(max_read_cov/10)) {
    min_read_cov <- 0
  }
  startingcoverages <- seq((min_read_cov+quarter_read_cov), max_read_cov, bottomtotop_read_cov)
  pattern <- c(rep(min_read_cov, 10000/windowsize), rep(startingcoverages[1], nrow(viral_subset)-(20000/windowsize)),rep(min_read_cov, 10000/windowsize))
  diff <- mean(abs(Cov_values_contig - pattern))
  start_pos <- (which(pattern == max(pattern))[1])
  end_pos <- which(pattern==max(pattern))[length(which(pattern==max(pattern)))]
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], "NA", start_pos, end_pos)
  for(cov in startingcoverages) {
    pattern <- c(rep(min_read_cov, 10000/windowsize), rep(cov, nrow(viral_subset)-(20000/windowsize)), rep(min_read_cov, 10000/windowsize))
    repeat {
      middle_rows <- which(pattern == cov)
      if (length(middle_rows) < blocksize/windowsize) break
      best_match_info <- blockpattern_translator(viral_subset, best_match_info, windowsize, pattern)
      pattern <- c(pattern[-c(middle_rows[2]:middle_rows[(2000/windowsize)+1])],rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
    }
  }
  new_pattern_coverages <- blockheight_optim(best_match_info, startingcoverages, bottomtotop_read_cov)
  for(newcov in new_pattern_coverages) {
    pattern <- c(rep(min_read_cov, 10000/windowsize), rep(newcov, nrow(viral_subset)-(20000/windowsize)), rep(min_read_cov, 10000/windowsize))
    repeat {
      best_match_info <- blockpattern_translator(viral_subset, best_match_info, windowsize, pattern)
      middle_rows <- which(pattern == newcov)
      pattern <- c(pattern[-c(middle_rows[2]:middle_rows[(2000/windowsize)+1])],rep(min_read_cov,2000/windowsize)) #remove 2000bp at a time
      if (length(middle_rows) < blocksize/windowsize) break #variable,
    }
  }
  best_match_results <- append(best_match_info, "Prophage-like")
  return(best_match_results)
}
