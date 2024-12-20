#' Sloping pattern builder
#'
#' Build a sloping pattern that consists of a sloping line spanning the contig
#' being assessed. The line slopes from left to right. The slope of the line is
#' changed, but the pattern is not translated across the contig.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileup
#' @param minSlope The minimum slope value to test for sloping patterns
#' @return List containing two objects
#' @keywords internal
fullSlope <- function(viralSubset, windowSize, minSlope) {
  maxReadCov <- max(viralSubset[, 2])
  minReadCov <- min(viralSubset[, 2])
  halfReadCov <- abs(maxReadCov - minReadCov) / 2
  newmax <-
    maxReadCov + ((abs(maxReadCov - (minReadCov + halfReadCov)) / 10))
  halfToMaxReadCov <-
    abs((newmax - (minReadCov + halfReadCov)) / 10)
  bestMatchInfoLR <-
    makeFullSlopes("Left", viralSubset, newmax, minReadCov, windowSize)
  bestMatchInfoRL <-
    makeFullSlopes("Right", viralSubset, newmax, minReadCov, windowSize)
  lapply(
    seq(newmax, (minReadCov + halfReadCov), -halfToMaxReadCov),
    function(cov) {
      slopeBottom <- minReadCov
      repeat {
        slopeChangeInfoL <-
          changeSlope(
            "Left",
            slopeBottom,
            halfToMaxReadCov,
            cov,
            viralSubset,
            windowSize
          )
        slopeChangeInfoR <-
          changeSlope(
            "Right",
            slopeBottom,
            halfToMaxReadCov,
            cov,
            viralSubset,
            windowSize
          )
        if (abs(slopeChangeInfoL[[4]] / windowSize) <
          minSlope | slopeChangeInfoL[[4]] / windowSize > 0) {
          break
        }
        if (abs(slopeChangeInfoR[[4]] / windowSize) <
          minSlope | slopeChangeInfoR[[4]] / windowSize < 0) {
          break
        }
        if (slopeChangeInfoL[[1]] < bestMatchInfoLR[[1]]) {
          bestMatchInfoLR <<- slopeChangeInfoL
        }
        if (slopeChangeInfoR[[1]] < bestMatchInfoRL[[1]]) {
          bestMatchInfoRL <<- slopeChangeInfoR
        }
        slopeBottom <- slopeChangeInfoL[[2]]
        slopeBottom <- slopeChangeInfoL[[2]]
      }
    }
  )
  return(list(bestMatchInfoLR, bestMatchInfoRL))
}


#' Make full slope patterns
#'
#' Makes slope patterns sloping either left to right (Left) or right to left
#' (right) across the contig being assessed.
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e.
#'   'Left') or positive slope (right to left, i.e. 'Right')
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param newMax A value for the top of the sloping pattern that is slightly
#'   higher than the maximum coverage value on the viralSubset
#' @param minReadCov Minimum read coverage value of the viralSubset
#' @param windowSize The window size used to re-average read coverage pileups
#' @return List
#' @keywords internal
makeFullSlopes <-
  function(leftOrRight,
           viralSubset,
           newMax,
           minReadCov,
           windowSize) {
    contigCoverage <- viralSubset[, 2]
    covSteps <- (newMax - minReadCov) / ((nrow(viralSubset) - 1))
    covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
    startSeq <- ifelse(leftOrRight == "Left", newMax, minReadCov)
    endSeq <- ifelse(leftOrRight == "Left", minReadCov, newMax)
    pattern <- seq(startSeq, endSeq, covSteps)
    diff <- mean(abs(contigCoverage - pattern))
    return(list(
      diff,
      minReadCov,
      newMax,
      covSteps,
      1,
      length(pattern),
      "Sloping"
    ))
  }

#' Change slope of sloping pattern
#'
#' Change the value of the slope used for the sloping pattern-match
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e.
#'   'Left') or positive slope (right to left, i.e. 'Right')
#' @param slopeBottom The value for the bottom of the sloping value
#' @param halfToMaxReadCov Half of the max VLP-fraction read coverage divided by
#'   10
#' @param cov The value for the top of the slope
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileup
#' @return List
#' @keywords internal
changeSlope <-
  function(leftOrRight,
           slopeBottom,
           halfToMaxReadCov,
           cov,
           viralSubset,
           windowSize) {
    slopeBottom <- slopeBottom + halfToMaxReadCov
    covSteps <- (cov - slopeBottom) / ((nrow(viralSubset) - 1))
    covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
    pattern <-
      if (leftOrRight == "Left") {
        (seq(cov, slopeBottom, covSteps))
      } else {
        (seq(slopeBottom, cov, covSteps))
      }
    diff <- mean(abs(viralSubset[, 2] - pattern))
    return(list(
      diff,
      slopeBottom,
      cov,
      covSteps,
      1,
      length(pattern),
      "Sloping"
    ))
  }
