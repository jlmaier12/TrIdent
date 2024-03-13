#' Build sloping pattern that slopes left-to-right with an initial jump-up in read coverage
#'
#' Build a sloping pattern that consists of a sloping line that starts on the contig being assessed. The pattern is used as input for the slopepattern_translator. The slope of the pattern is changed after each full translation across a contig.
#'
#' @param viral_subset A subset of the read coverage dataset that pertains only to the contig currently being assessed
#' @param windowsize The window size used to re-average read coverage datasets
#' @keywords internal
slope_LefttoRight_withstart <- function (viral_subset, windowsize) {
  max_read_cov <- max(viral_subset[,2])
  min_read_cov <- min(viral_subset[,2])
  half_read_cov <- abs((max_read_cov-min_read_cov))/2
  bottomtotop_read_cov <- as.numeric(abs((max_read_cov - (min_read_cov+half_read_cov)))/10)
  newmax <- max_read_cov+bottomtotop_read_cov
  Cov_values_contig <- viral_subset[,2]
  cov_steps <- (newmax-min_read_cov)/((nrow(viral_subset)-((10000/windowsize)+1))) 
  pattern <- c(rep(min_read_cov,10000/windowsize), seq(newmax,min_read_cov, -cov_steps)) #different in each function
  diff <- mean(abs(Cov_values_contig - pattern))
  start_pos <- which(pattern==max(pattern)) #different in each function
  slope <- (min_read_cov-newmax)/(nrow(viral_subset)-1) #different in each function
  best_match_info <- list(diff, min_read_cov, newmax, -cov_steps, start_pos, length(pattern), slope) #different in each function
  lapply(seq(newmax, (min_read_cov+half_read_cov), -bottomtotop_read_cov), function(cov) {
    slope_bottom <- min_read_cov
    cov_steps <- (cov-slope_bottom)/((nrow(viral_subset)-((10000/windowsize)+1)))
    pattern <- c(rep(min_read_cov,10000/windowsize), seq(cov,slope_bottom, -cov_steps)) #different in each function
    slope <- (cov-slope_bottom)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
    step <- ((slope_bottom-cov)/10) #different in each function
    if (abs(slope) < (15/100000) | slope >0) return(NULL) #different in each function
    repeat {
      best_match_info <- slopepattern_translator(viral_subset, best_match_info, windowsize, pattern, "lefttoright") #different in each function
      slope_bottom <- slope_bottom + step
      cov_steps <- (cov-slope_bottom)/((nrow(viral_subset)-((10000/windowsize)+1)))
      pattern <- c(rep(min_read_cov,10000/windowsize), seq(cov,slope_bottom, -cov_steps)) #different in each function
      slope <- (slope_bottom-cov)/(nrow(viral_subset)-((10000/windowsize)+1)) #different in each function
      if (abs(slope) < (15/100000) | slope >0) break #different in each function
    }
  })
  best_match_results <- c(best_match_info, "Gen/Lat/GTA")
  return(best_match_results)
}
