#' Make slope patterns with starts
#'
#' Makes slope patterns sloping either left to right (Left) or right to left (right) across the contig being assessed.
#' Slope patterns contain an initiation point.
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e. 'Left') or positive slope (right to left, i.e. 'Right')
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param newMax A value for the top of the sloping pattern that is slightly higher than the maximum coverage value on the viralSubset
#' @param minReadCov Minimum read coverage value of the viralSubset
#' @param windowSize The window size used to re-average read coverage pileup
#' @keywords internal
makeSlopesWStarts <- function(leftOrRight, viralSubset, newMax, minReadCov, windowSize){
  contigCoverage <- viralSubset[,2]
  covSteps <- (newMax - minReadCov) / ((nrow(viralSubset) - ((5000 / windowSize) + 1)))
  covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
  pattern <- if(leftOrRight == "Left") c(rep(minReadCov, 5000 / windowSize), seq(newMax, minReadCov, covSteps))
             else c(seq(minReadCov, newMax, covSteps), rep(minReadCov, 5000 / windowSize))
  diff <- mean(abs(contigCoverage - pattern))
  slope <-  ifelse(leftOrRight == "Left", (minReadCov - newMax) / ((nrow(viralSubset) * windowSize) - 5000), (newMax - minReadCov) /  ((nrow(viralSubset) * windowSize) - 5000))
  startPos <- ifelse(leftOrRight == "Left", which(pattern == max(pattern)), 1)
  endPos <- ifelse(leftOrRight == "Left", length(pattern), which(pattern == max(pattern)))
  return(list(diff, minReadCov, newMax, covSteps, startPos, endPos, slope, "Sloping"))
}
