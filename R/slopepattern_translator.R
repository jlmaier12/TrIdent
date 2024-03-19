#' Sloping pattern translator
#'
#' Translates a sloping pattern containing the initial jump-up in read coverage across a contig. Translate the pattern 2000bp at a time. Stop translating when the pattern left on the contig reaches 45,000bp.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param best_match_info The information associated with the current best pattern match. Includes the match score, the min and max y-axis values, the length and the start and stop positions of the pattern
#' @param windowsize The window size used to re-average read coverage datasets
#' @param pattern A vector containing the values associated with the generalized pattern to be translated across the contig
#' @param direction The direction of the slope of the generalized pattern being translated across. Either "lefttoright" or "righttoleft"
#' @keywords internal
slopepattern_translator <- function(viral_subset, best_match_info, windowsize, pattern, direction){
  min_read_cov <- min(pattern)
  repeat {
    if (direction== "lefttoright") {
      pattern <- c(rep(min_read_cov,2000/windowsize),pattern[-c((length(pattern)-((2000/windowsize)-1)):length(pattern))])
      start_pos <- which(pattern==max(pattern))
      slope_bottom <- min(pattern[pattern != min(pattern)])
      end_pos <- which(pattern==slope_bottom)
      cov_steps <- -(max(pattern)-slope_bottom)/abs(end_pos-start_pos)
    }
    if (direction =="righttoleft") {
      pattern <- c(pattern[-c(1:(2000/windowsize))],rep(min_read_cov,2000/windowsize))
      end_pos <- which(pattern==max(pattern))
      slope_bottom <- min(pattern[pattern != min(pattern)])
      start_pos <- which(pattern==slope_bottom)
      cov_steps <- ((max(pattern)-slope_bottom)/abs(end_pos-start_pos))
    }
    if((length(pattern[!(pattern %in% min_read_cov)]) * windowsize)<45000) break
    diff <- mean(abs(viral_subset[,2] - pattern))
    if (diff < best_match_info[[1]]){
      best_match_info <- list(diff, slope_bottom, max(pattern), cov_steps, start_pos, end_pos, best_match_info[[7]])
    }
  }
  return(best_match_info)
}
