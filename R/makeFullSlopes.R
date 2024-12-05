#' Make full slope patterns
#'
#' Makes slope patterns sloping either left to right (Left) or right to left (right) across the contig being assessed.
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e. 'Left') or positive slope (right to left, i.e. 'Right')
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param newMax A value for the top of the sloping pattern that is slightly higher than the maximum coverage value on the viralSubset
#' @param minReadCov Minimum read coverage value of the viralSubset
#' @param windowSize The window size used to re-average read coverage pileups
#' @keywords internal
makeFullSlopes <- function(leftOrRight, viralSubset, newMax, minReadCov, windowSize){
  contigCoverage <- viralSubset[,2]
  covSteps <- (newMax - minReadCov) / ((nrow(viralSubset) - 1))
  covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
  startSeq <- ifelse(leftOrRight == "Left", newMax, minReadCov)
  endSeq <- ifelse(leftOrRight == "Left", minReadCov, newMax)
  pattern <- seq(startSeq, endSeq, covSteps)
  diff <- mean(abs(contigCoverage - pattern))
  slope <- (endSeq - startSeq) / (nrow(viralSubset) * windowSize)
  return(list(diff, minReadCov, newMax, covSteps, 1, length(pattern), slope, "Sloping"))
}
