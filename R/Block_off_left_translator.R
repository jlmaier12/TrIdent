#' Block off-left pattern translator
#'
#' Build and translate a block pattern going off the left side of the contig. Stop translating when the pattern left on the contig is 10,000bp. Translate the pattern 2000bp at a time.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param blocksize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @keywords internal
block_off_left_translator <- function (viral_subset, windowsize, blocksize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  quarter_read_cov <- abs(max_read_cov-min_read_cov)/4
  bottomtotop_read_cov <- abs(max_read_cov - (min_read_cov+quarter_read_cov))/10
  Cov_values_contig <- viral_subset[,2]
  if(min_read_cov>(max_read_cov/10)) {
    min_read_cov <- 0
  }
  startingcoverages <- seq((min_read_cov+quarter_read_cov),  max_read_cov, bottomtotop_read_cov)
  pattern <- c(rep(startingcoverages[1], nrow(viral_subset)-(10000/windowsize)), rep(min_read_cov, 10000/windowsize))
  diff <- mean(abs(Cov_values_contig - pattern))
  end_pos <- (which(pattern == min(pattern))[1])-1
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], "NA", 1, end_pos)
  for(cov in startingcoverages) {
    pattern <- c(rep(cov, nrow(viral_subset)-(10000/windowsize)), rep(min_read_cov, 10000/windowsize))
    repeat {
      diff <- mean(abs(Cov_values_contig - pattern))
      end_pos <- (which(pattern == min_read_cov)[1])-1
      if (diff < best_match_info[[1]]){
        best_match_info <- list(diff, min_read_cov, cov, "NA", 1, end_pos)
      }
      pattern <- c(pattern[-c(1:(2000/windowsize))], rep(min_read_cov, (2000/windowsize)))
      if (length(which(pattern==cov)) < blocksize/windowsize) break
    }
  }
  new_pattern_coverages <- blockheight_optim(best_match_info, startingcoverages, bottomtotop_read_cov)
  for(newcov in new_pattern_coverages) {
    pattern <- c(rep(newcov,  nrow(viral_subset)-(10000/windowsize)), rep(min_read_cov, 10000/windowsize))
    repeat {
      end_pos <- (which(pattern == min(pattern))[1])-1
      diff <- mean(abs(Cov_values_contig - pattern))
      if (diff < best_match_info[[1]]){
        best_match_info <-  list(diff, min_read_cov, cov, "NA", 1, end_pos)
      }
      pattern <- c(pattern[-c(1:(2000/windowsize))], rep(min_read_cov, (2000/windowsize)))
      if (length(which(pattern==newcov)) < blocksize/windowsize) break
    }
  }
  best_match_results <- append(best_match_info, "Prophage-like")
  return(best_match_results)
}

