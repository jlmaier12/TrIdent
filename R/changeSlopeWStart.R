#' Change slope of sloping pattern with initial start
#'
#' Change the value of the slope used for the sloping with start pattern-match
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e. 'Left') or positive slope (right to left, i.e. 'Right')
#' @param slopeBottom The value for the bottom of the sloping value
#' @param slopeBottomChange The value used to increase the bottom of the slope
#' @param cov The value for the top of the slope
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileup
#' @keywords internal
#'
changeSlopeWStart <- function(leftOrRight, slopeBottom, slopeBottomChange, cov, viralSubset, windowSize){
  minReadCov <- min(viralSubset[,2])
  slopeBottom <- slopeBottom + slopeBottomChange
  covSteps <- (cov - slopeBottom) / ((nrow(viralSubset) - ((5000 / windowSize) + 1)))
  covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
  pattern <- if(leftOrRight == "Left") c(rep(minReadCov, 5000 / windowSize), seq(cov, slopeBottom, covSteps))
             else c(seq(slopeBottom, cov, covSteps), rep(minReadCov, 5000 / windowSize))
  slope <-  ifelse(leftOrRight == "Left", (slopeBottom - cov) / ((nrow(viralSubset) * windowSize) - 5000), (cov - slopeBottom) /  ((nrow(viralSubset) * windowSize) - 5000))
  return(list(pattern, slope, slopeBottom))
}

