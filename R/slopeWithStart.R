#' Sloping pattern with an initial jump-up in read coverage
#'
#' Build, translate, and change slope of sloping pattern with slope start
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileups
#' @keywords internal
slopeWithStart <- function (viralSubset, windowSize) {
  maxReadCov <- max(viralSubset[,2])
  minReadCov <- min(viralSubset[,2])
  halfReadCov <- abs((maxReadCov - minReadCov)) / 2
  newMax <- maxReadCov + (abs((maxReadCov - (minReadCov + halfReadCov)) / 10))
  halfToMaxReadCov <- abs((newMax - (minReadCov + halfReadCov)) / 10)
  bestMatchInfoLR <- makeSlopesWStarts("Left", viralSubset, newMax, minReadCov, windowSize)
  bestMatchInfoRL <- makeSlopesWStarts("Right", viralSubset, newMax, minReadCov, windowSize)
  lapply(seq(newMax, (minReadCov + halfReadCov), -halfToMaxReadCov), function(cov) {
    slopeBottom <- minReadCov
    slopeBottomChange <- (cov - minReadCov) / 10
    repeat {
      slopeChangeLR <- changeSlopeWStart("Left", slopeBottom, slopeBottomChange, cov, viralSubset, windowSize)
      slopeChangeRL <- changeSlopeWStart("Right", slopeBottom, slopeBottomChange, cov, viralSubset, windowSize)
      if (abs(slopeChangeLR[[2]]) < (15 / 100000) | slopeChangeLR[[2]] > 0) break
      if (abs(slopeChangeRL[[2]]) < (15 / 100000) | slopeChangeRL[[2]] < 0) break
      bestMatchInfoLR <<- slopeTranslator(viralSubset, bestMatchInfoLR, windowSize, slopeChangeLR, "Left")
      bestMatchInfoRL <<- slopeTranslator(viralSubset, bestMatchInfoRL, windowSize, slopeChangeRL, "Right")
      slopeBottom <- slopeChangeLR[[3]]
    }
  })
  return(list(bestMatchInfoLR, bestMatchInfoRL))
}
