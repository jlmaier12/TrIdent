#' Block pattern height optimization
#'
#' Optimizes heights of block patterns based on scores of initial pattern matches
#'
#' @param best_match_info The information associated with the current best pattern match. Includes the match score, the min and max y-axis values, the length and the start and stop positions of the pattern
#' @param oldcoverages The initial coverage values used for pattern matching
#' @param oldcovstep The initial difference between coverage values used for pattern matching
#'
#' @keywords internal
blockheight_optim <- function (best_match_info, oldcoverages, oldcovstep){
  min_diff <- best_match_info[[1]]
  min_pattern_cov <- best_match_info[[2]]
  max_pattern_cov <- best_match_info[[3]]
  min_diff_index <- which(oldcoverages==max_pattern_cov)
  if (min_diff_index == 1) {
    new_cov_step <- abs(max_pattern_cov-(max_pattern_cov+oldcovstep))/5
    optimized_pattern_coverages <- seq(max_pattern_cov, (max_pattern_cov+oldcovstep), new_cov_step)
  } else if (min_diff_index==11) {
    new_cov_step <- abs((max_pattern_cov-oldcovstep)-max_pattern_cov)/5
    optimized_pattern_coverages <- seq((max_pattern_cov-oldcovstep), max_pattern_cov, new_cov_step)
  } else {
    new_cov_step <- abs((max_pattern_cov-oldcovstep)-(max_pattern_cov+oldcovstep))/10
    optimized_pattern_coverages <- seq((max_pattern_cov-oldcovstep), (max_pattern_cov+oldcovstep), new_cov_step)
  }
  return(optimized_pattern_coverages)
}
