#' Builds prophage-like block patterns
#'
#' Build and translate a block pattern going off the left side, right side and full length of the contig.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileups
#' @param minBlockSize The minimum size of the prophage-like block pattern. Default is 10000 bp.
#' @param maxBlockSize The maximum size of the prophage-like block pattern. Default is NA.
#' @return List containing three objects
#' @keywords internal
blockBuilder <- function (viralSubset, windowSize, minBlockSize, maxBlockSize) {
  maxReadCov <- max(viralSubset[,2])
  minReadCov <- min(viralSubset[,2])
  if(minReadCov > (maxReadCov * 0.2)) minReadCov <- 0
  quarterReadCov <- abs(maxReadCov - minReadCov) / 4
  covSteps <- abs(maxReadCov - (minReadCov + quarterReadCov)) / 10
  startingCovs <- seq((minReadCov + quarterReadCov),  maxReadCov, covSteps)

  blockLength <- ifelse((nrow(viralSubset) - (5000 / windowSize)) > (maxBlockSize / windowSize),
                        maxBlockSize / windowSize, nrow(viralSubset) - (5000 / windowSize))
  nonBlock <- nrow(viralSubset) - blockLength
  blockLengthFull <- ifelse((nrow(viralSubset) - (10000 / windowSize)) > (maxBlockSize / windowSize),
                            maxBlockSize / windowSize, nrow(viralSubset) - (10000 / windowSize))
  nonBlockFull <- nrow(viralSubset) - (blockLengthFull + (5000 / windowSize))

  bestMatchInfoFull <- makeBlockPattern(viralSubset, windowSize, "Full", blockLengthFull, nonBlockFull, minReadCov, startingCovs[1])[[1]]
  bestMatchInfoRight <- makeBlockPattern(viralSubset, windowSize, "Right", blockLength, nonBlock, minReadCov, startingCovs[1])[[1]]
  bestMatchInfoLeft <- makeBlockPattern(viralSubset, windowSize, "Left", blockLength, nonBlock, minReadCov, startingCovs[1])[[1]]

  lapply(seq_along(startingCovs), function(i) {
    cov <- startingCovs[[i]]
    patternFull <- makeBlockPattern(viralSubset, windowSize, "Full", blockLengthFull, nonBlockFull, minReadCov, cov)[[2]]
    patternRight <- makeBlockPattern(viralSubset, windowSize, "Right", blockLength, nonBlock, minReadCov, cov)[[2]]
    patternLeft <- makeBlockPattern(viralSubset, windowSize, "Left", blockLength, nonBlock, minReadCov, cov)[[2]]
    bestMatchInfoLeft <<- leftRightBlockTranslater(viralSubset, patternLeft, "Left", windowSize, minReadCov, cov, bestMatchInfoLeft, minBlockSize)
    bestMatchInfoRight <<- leftRightBlockTranslater(viralSubset, patternRight, "Right", windowSize, minReadCov, cov, bestMatchInfoRight, minBlockSize)
    repeat {
      bestMatchInfoFull <<- blockTranslator(viralSubset, bestMatchInfoFull, windowSize, patternFull)
      middleRows <- which(patternFull == cov)
      if (length(middleRows) < (minBlockSize / windowSize) + 1) break
      patternFull <- c(patternFull[-c(middleRows[2]:middleRows[(1000 / windowSize) + 1])], rep(minReadCov, 1000 / windowSize))
    }
  })
  return(list(bestMatchInfoLeft, bestMatchInfoRight, bestMatchInfoFull))
}
