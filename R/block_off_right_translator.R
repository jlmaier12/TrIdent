#' Block off-right pattern translator
#'
#' Build and translate a block pattern going off the right side of the contig. Stop translating when the pattern left on the contig is 10,000bp. Translate the pattern 2000bp at a time.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @param minblocksize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @param maxblocksize The maximum size of the prophage-like block pattern. Default is NA.
#' @keywords internal
block_off_right_translator <- function (viral_subset, windowsize, minblocksize, maxblocksize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  quarter_read_cov <- abs(max_read_cov-min_read_cov)/4
  bottomtotop_read_cov <- abs(max_read_cov - (min_read_cov+quarter_read_cov))/10
  Cov_values_contig <- viral_subset[,2]
  if(min_read_cov>(max_read_cov*0.2)) {
    min_read_cov <- 0
  }
  startingcoverages <- seq((min_read_cov+quarter_read_cov),  max_read_cov, bottomtotop_read_cov)
  shape_length <- ifelse ((nrow(viral_subset)-(10000/windowsize))>(maxblocksize/windowsize),maxblocksize/windowsize,nrow(viral_subset)-(10000/windowsize))
  nonshape <- nrow(viral_subset)-shape_length
  pattern <- c(rep(min_read_cov, nonshape), rep(startingcoverages[1], shape_length)) #Different in each function
  diff <- mean(abs(Cov_values_contig - pattern))
  start_pos <- (which(pattern == max(pattern))[1]) #The start and stop positions used for determining the pattern length are different in each function
  best_match_info <- list(diff, min_read_cov, startingcoverages[1], "NA", start_pos, length(pattern), "NA") #Recorded different for each function
  lapply(1:length(startingcoverages), function(i) {
    cov<-startingcoverages[[i]]
    pattern <- c(rep(min_read_cov, nonshape), rep(cov, shape_length)) #Different in each function
    repeat {
      diff <- mean(abs(Cov_values_contig - pattern))
      start_pos <- (which(pattern == max(pattern))[1])
      if (diff < best_match_info[[1]]){
        best_match_info <<- list(diff, min_read_cov, cov, "NA", start_pos, length(pattern), "NA") #Different in each function
      }
      pattern <- c(rep(min_read_cov,(2000/windowsize)),pattern[-c(((length(pattern))-((2000/windowsize)-1)):length(pattern))]) #Different in each function
      if (length(which(pattern==cov)) < (minblocksize/windowsize)+1) break
    }
  })
  new_pattern_coverages <- blockheight_optim(best_match_info, startingcoverages, bottomtotop_read_cov)
  lapply(1:length(new_pattern_coverages), function(i) {
    newcov<-new_pattern_coverages[[i]]
    pattern <- c(rep(min_read_cov, nonshape), rep(newcov, shape_length)) #Different in each function
    repeat {
      diff <- mean(abs(Cov_values_contig - pattern))
      start_pos <- (which(pattern == max(pattern))[1]) #Different in each function
      if (diff < best_match_info[[1]]){
        best_match_info <<- list(diff, min_read_cov, newcov, "NA", start_pos, length(pattern), "NA")
      }
      pattern <- c(rep(min_read_cov,(2000/windowsize)),pattern[-c(((length(pattern))-((2000/windowsize)-1)):length(pattern))]) ##Different in each function
      if (length(which(pattern==newcov)) < (minblocksize/windowsize)+1) break
    }
  })
  best_match_results <- c(best_match_info, "Prophage-like")
  return(best_match_results)
}
