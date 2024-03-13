#' Sloping pattern right-to-left pattern builder
#'
#' Build a sloping pattern that consists of a sloping line spanning the contig being assessed. The line slopes from right to left. The slope of the line is changed, but the pattern is not translated across the contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize <- The window size used to re-average read coverage datasets
#' @keywords internal
slope_RighttoLeft <- function (viral_subset, windowsize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  half_read_cov <- abs((max_read_cov-min_read_cov))/2
  bottomtotop_read_cov <- as.numeric(abs((max_read_cov - (min_read_cov+half_read_cov)))/10)
  newmax <- max_read_cov+bottomtotop_read_cov
  bottomtotop_read_cov <- abs((newmax-(min_read_cov+half_read_cov))/10)
  Cov_values_contig <- viral_subset[,2]
  cov_steps <- (newmax-min_read_cov)/((nrow(viral_subset)-1))
  pattern <- seq(min_read_cov,newmax,cov_steps) #different in each function
  diff <- mean(abs(Cov_values_contig - pattern))
  slope <- (newmax-min_read_cov)/(nrow(viral_subset)-1) #different in each function
  best_match_info <- list(diff, min_read_cov, newmax, cov_steps, 1, length(pattern), slope) #different in each function
  lapply(seq(newmax, (min_read_cov+half_read_cov), -bottomtotop_read_cov), function(cov) {
    slope_bottom <- min_read_cov
    cov_steps <- (cov-slope_bottom)/((nrow(viral_subset)-1))
    pattern <- seq(slope_bottom,cov,cov_steps) #different in each function
    slope <- (cov-slope_bottom)/(nrow(viral_subset)-1) #different in each function
    if (abs(slope) < 15/100000 | slope < 0) return(NULL)
    repeat {
      if (diff < best_match_info[[1]]) {
        best_match_info <<- list(diff, slope_bottom, cov, cov_steps, 1, length(pattern), slope) #different in each function
      }
      slope_bottom <- slope_bottom + bottomtotop_read_cov
      cov_steps <- (cov-slope_bottom)/((nrow(viral_subset)-1))
      pattern <- seq(slope_bottom,cov,cov_steps) #different in each function
      diff <- mean(abs(Cov_values_contig - pattern))
      slope <- (cov-slope_bottom)/(nrow(viral_subset)-1) #different in each function
      if (abs(slope) < 15/100000 | slope < 0) break #different in each function
    }
  })
  best_match_results <- c(best_match_info, "Gen/Lat/GTA")
  return(best_match_results)
}
