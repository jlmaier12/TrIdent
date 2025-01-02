#' Sloping pattern with an initial jump-up in read coverage
#'
#' Build, translate, and change slope of sloping pattern with slope start
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileups
#' @param minSlope The minimum slope value to test for sloping patterns
#' @return List containing two objects
#' @keywords internal
slopeWithStart <- function(viralSubset, windowSize, minSlope) {
  maxReadCov <- max(viralSubset[, 2])
  minReadCov <- min(viralSubset[, 2])
  halfReadCov <- abs((maxReadCov - minReadCov)) / 2
  newMax <-
    maxReadCov + (abs((maxReadCov - (
      minReadCov + halfReadCov
    )) / 10))
  halfToMaxReadCov <-
    abs((newMax - (minReadCov + halfReadCov)) / 10)
  bestMatchInfoLR <- makeSlopesWStarts(
    "Left", viralSubset, newMax,
    minReadCov, windowSize
  )
  bestMatchInfoRL <- makeSlopesWStarts(
    "Right", viralSubset, newMax,
    minReadCov, windowSize
  )
  for(cov in seq(newMax, (minReadCov + halfReadCov), -halfToMaxReadCov)){
      slopeBottom <- minReadCov
      slopeBottomChange <- (cov - minReadCov) / 10
      repeat {
        slopeChangeLR <- changeSlopeWStart(
          "Left",
          slopeBottom,
          slopeBottomChange,
          cov,
          viralSubset,
          windowSize
        )
        slopeChangeRL <- changeSlopeWStart(
          "Right",
          slopeBottom,
          slopeBottomChange,
          cov,
          viralSubset,
          windowSize
        )
        if (abs(slopeChangeLR[[2]] / windowSize) < minSlope |
          slopeChangeLR[[2]] / windowSize > 0) {
          break
        }
        if (abs(slopeChangeRL[[2]] / windowSize) < minSlope |
          slopeChangeRL[[2]] / windowSize < 0) {
          break
        }
        bestMatchInfoLR <-
          slopeTranslator(
            viralSubset,
            bestMatchInfoLR,
            windowSize,
            slopeChangeLR,
            "Left"
          )
        bestMatchInfoRL <-
          slopeTranslator(
            viralSubset,
            bestMatchInfoRL,
            windowSize,
            slopeChangeRL,
            "Right"
          )
        slopeBottom <- slopeChangeLR[[3]]
      }
    }
  return(list(bestMatchInfoLR, bestMatchInfoRL))
}


#' Make slope patterns with starts
#'
#' Makes slope patterns sloping either left to right (Left) or right to left
#' (right) across the contig being assessed. Slope patterns contain an
#' initiation point.
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e.
#'   'Left') or positive slope (right to left, i.e. 'Right')
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param newMax A value for the top of the sloping pattern that is slightly
#'   higher than the maximum coverage value on the viralSubset
#' @param minReadCov Minimum read coverage value of the viralSubset
#' @param windowSize The window size used to re-average read coverage pileup
#' @return List
#' @keywords internal
makeSlopesWStarts <-
  function(leftOrRight,
           viralSubset,
           newMax,
           minReadCov,
           windowSize) {
    contigCoverage <- viralSubset[, 2]
    covSteps <- (newMax - minReadCov) / ((nrow(viralSubset) -
      ((5000 / windowSize) + 1)))
    covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
    pattern <-
      if (leftOrRight == "Left") {
        c(
          rep(minReadCov, 5000 / windowSize),
          seq(newMax, minReadCov, covSteps)
        )
      } else {
        c(
          seq(minReadCov, newMax, covSteps),
          rep(minReadCov, 5000 / windowSize)
        )
      }
    diff <- mean(abs(contigCoverage - pattern))
    startPos <- ifelse(leftOrRight == "Left",
      which(pattern == max(pattern)), 1
    )
    endPos <- ifelse(leftOrRight == "Left",
      length(pattern), which(pattern == max(pattern))
    )
    return(list(
      diff,
      minReadCov,
      newMax,
      covSteps,
      startPos,
      endPos,
      "Sloping"
    ))
  }

#' Change slope of sloping pattern with initial start
#'
#' Change the value of the slope used for the sloping with start pattern-match
#'
#' @param leftOrRight Generate pattern for negative slope (left to right, i.e.
#'   'Left') or positive slope (right to left, i.e. 'Right')
#' @param slopeBottom The value for the bottom of the sloping value
#' @param slopeBottomChange The value used to increase the bottom of the slope
#' @param cov The value for the top of the slope
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileup
#' @return List
#' @keywords internal
#'
changeSlopeWStart <-
  function(leftOrRight,
           slopeBottom,
           slopeBottomChange,
           cov,
           viralSubset,
           windowSize) {
    minReadCov <- min(viralSubset[, 2])
    slopeBottom <- slopeBottom + slopeBottomChange
    covSteps <-
      (cov - slopeBottom) / ((nrow(viralSubset) - ((5000 / windowSize)
      + 1)))
    covSteps <- ifelse(leftOrRight == "Left", -covSteps, covSteps)
    pattern <-
      if (leftOrRight == "Left") {
        c(
          rep(minReadCov, 5000 / windowSize),
          seq(cov, slopeBottom, covSteps)
        )
      } else {
        c(
          seq(slopeBottom, cov, covSteps),
          rep(minReadCov, 5000 / windowSize)
        )
      }
    return(list(pattern, covSteps, slopeBottom))
  }



#' Sloping pattern translator
#'
#' Translates a sloping pattern containing the initial jump-up in read coverage
#' across a contig. Translate the pattern 1000 bp at a time. Stop translating
#' when the pattern left on the contig reaches 20,000 bp.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param bestMatchInfo The pattern-match information associated with the
#'   current best pattern match.
#' @param windowSize The window size used to re-average read coverage pileups
#' @param slopeChange A list containing pattern vector, slope value, and
#'   value of slope bottom
#' @param leftOrRight The direction of the sloping pattern. Either "Left" for
#'   left to right (neg) slopes or "Right" for right to left (pos) slopes.
#' @return List
#' @keywords internal
slopeTranslator <-
  function(viralSubset,
           bestMatchInfo,
           windowSize,
           slopeChange,
           leftOrRight) {
    pattern <- slopeChange[[1]]
    minPattern <- min(pattern)
    repeat {
      if (leftOrRight == "Left") {
        pattern <- c(
          rep(minPattern, 1000 / windowSize),
          pattern[-c((length(pattern) - ((
            1000 / windowSize
          ) - 1)):length(pattern))]
        )
        slopeBottomIdx <-
          min(pattern[pattern != min(pattern)])
        startRowIdx <- which(pattern == max(pattern))
        endRowIdx <- which(pattern == slopeBottomIdx)
      }
      if (leftOrRight == "Right") {
        pattern <- c(
          pattern[-(seq_len(1000 / windowSize))],
          rep(minPattern, 1000 / windowSize)
        )
        slopeBottomIdx <-
          min(pattern[pattern != min(pattern)])
        startRowIdx <- which(pattern == slopeBottomIdx)
        endRowIdx <- which(pattern == max(pattern))
      }
      if ((length(pattern[!(pattern %in% minPattern)]) * windowSize) <
        20000) {
        break
      }
      diff <- mean(abs(viralSubset[, 2] - pattern))
      if (diff < bestMatchInfo[[1]]) {
        covSteps <- ((max(pattern) - slopeBottomIdx) /
          abs(endRowIdx - startRowIdx))
        covSteps <-
          ifelse(leftOrRight == "Left", -covSteps, covSteps)
        bestMatchInfo <-
          list(
            diff,
            slopeBottomIdx,
            max(pattern),
            covSteps,
            startRowIdx,
            endRowIdx,
            "Sloping"
          )
      }
    }
    return(bestMatchInfo)
  }
