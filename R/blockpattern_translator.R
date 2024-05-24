#' Full block pattern-translator
#'
#' Translates full block-pattern across a contig. Translate the pattern 2000bp at a time. Stop translating when the pattern is 5000bp from the end of the contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param best_match_info The information associated with the current best pattern match. Includes the match score, the min and max y-axis values, the length and the start and stop positions of the pattern
#' @param windowsize The window size used to re-average read coverage datasets
#' @param pattern A vector containing the values associated with the block pattern to be translated across the contig
#'
#' @keywords internal
blockpattern_translator <- function(viral_subset, best_match_info, windowsize, pattern){
  Cov_values_contig <- viral_subset[,2]
  min_pattern_cov <- min(pattern)
  max_pattern_cov <- max(pattern)
  repeat {
    pattern <- c(rep(min_pattern_cov, (1000/windowsize)),pattern[-c((length(pattern)-((1000/windowsize)-1)):length(pattern))])
    if(pattern[length(pattern)-(5000/windowsize)]>min_pattern_cov) break
    diff <- mean(abs(Cov_values_contig - pattern))
    start_pos <- which(pattern==max(pattern))[1]
    end_pos <- which(pattern==max(pattern))[length(which(pattern==max(pattern)))]
    if (diff < best_match_info[[1]]){
      best_match_info <- list(diff, min_pattern_cov, max_pattern_cov, "NA", start_pos, end_pos, "NA")
    }
  }
  return(best_match_info)
}
