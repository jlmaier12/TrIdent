#' Builds prophage-like block patterns
#'
#' Build and translate a block pattern going off the left side, right side and
#' full length of the contig.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileups
#' @param minBlockSize The minimum size of the prophage-like block pattern.
#'   Default is 10000 bp.
#' @param maxBlockSize The maximum size of the prophage-like block pattern.
#'   Default is NA.
#' @return List containing three objects
#' @keywords internal
blockBuilder <-
  function(viralSubset,
           windowSize,
           minBlockSize,
           maxBlockSize) {
    maxReadCov <- max(viralSubset[, 2])
    minReadCov <- min(viralSubset[, 2])
    if (minReadCov > (maxReadCov * 0.2)) {
      minReadCov <- 0
    }
    quarterReadCov <- abs(maxReadCov - minReadCov) / 4
    covSteps <- abs(maxReadCov - (minReadCov + quarterReadCov)) / 10
    startingCovs <-
      seq((minReadCov + quarterReadCov), maxReadCov, covSteps)
    blockLength <- ifelse((nrow(viralSubset) - (5000 / windowSize)) >
      (maxBlockSize / windowSize),
    maxBlockSize / windowSize,
    nrow(viralSubset) - (5000 / windowSize)
    )
    nonBlock <- nrow(viralSubset) - blockLength
    blockLengthFull <-
      ifelse((nrow(viralSubset) - (10000 / windowSize)) >
        (maxBlockSize / windowSize),
      maxBlockSize / windowSize,
      nrow(viralSubset) - (10000 / windowSize)
      )
    nonBlockFull <-
      nrow(viralSubset) - (blockLengthFull + (5000 / windowSize))
    bestMatchInfoFull <-
      makeBlockPattern(
        viralSubset,
        windowSize,
        "Full",
        blockLengthFull,
        nonBlockFull,
        minReadCov,
        startingCovs[1]
      )[[1]]
    bestMatchInfoRight <-
      makeBlockPattern(
        viralSubset,
        windowSize,
        "Right",
        blockLength,
        nonBlock,
        minReadCov,
        startingCovs[1]
      )[[1]]
    bestMatchInfoLeft <-
      makeBlockPattern(
        viralSubset,
        windowSize,
        "Left",
        blockLength,
        nonBlock,
        minReadCov,
        startingCovs[1]
      )[[1]]
    for (i in seq_along(startingCovs)) {
      cov <- startingCovs[[i]]
      patternFull <- makeBlockPattern(
        viralSubset,
        windowSize,
        "Full",
        blockLengthFull,
        nonBlockFull,
        minReadCov,
        cov
      )[[2]]
      patternRight <-
        makeBlockPattern(
          viralSubset,
          windowSize,
          "Right",
          blockLength,
          nonBlock,
          minReadCov,
          cov
        )[[2]]
      patternLeft <- makeBlockPattern(
        viralSubset,
        windowSize,
        "Left",
        blockLength,
        nonBlock,
        minReadCov,
        cov
      )[[2]]
      bestMatchInfoLeft <-
        leftRightBlockTranslater(
          viralSubset,
          patternLeft,
          "Left",
          windowSize,
          minReadCov,
          cov,
          bestMatchInfoLeft,
          minBlockSize
        )
      bestMatchInfoRight <-
        leftRightBlockTranslater(
          viralSubset,
          patternRight,
          "Right",
          windowSize,
          minReadCov,
          cov,
          bestMatchInfoRight,
          minBlockSize
        )
      repeat {
        bestMatchInfoFull <-
          blockTranslator(
            viralSubset, bestMatchInfoFull,
            windowSize, patternFull
          )
        middleRows <- which(patternFull == cov)
        if (length(middleRows) < (minBlockSize / windowSize) + 1) {
          break
        }
        patternFull <- c(
          patternFull[-c(middleRows[2]:
          middleRows[(1000 / windowSize) + 1])],
          rep(minReadCov, 1000 / windowSize)
        )
      }
    }
    return(list(bestMatchInfoLeft, bestMatchInfoRight, bestMatchInfoFull))
  }

#' Make block patterns for pattern-matching
#'
#' Make full, left and right block patterns for Prophage-like classifications
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param windowSize The window size used to re-average read coverage pileups
#' @param fullLeftRight The block pattern variation being built
#' @param blockLength Maximum block pattern length
#' @param nonBlock Maximum non-block pattern length
#' @param minReadCov Either 0 or the minimum VLP-fraction read coverage value
#' @param cov The height value of the block pattern
#' @return List containing two objects
#' @keywords internal
makeBlockPattern <- function(viralSubset,
                             windowSize,
                             fullLeftRight,
                             blockLength,
                             nonBlock,
                             minReadCov,
                             cov) {
  if (fullLeftRight == "Full") {
    pattern <- c(
      rep(minReadCov, 5000 / windowSize),
      rep(cov, blockLength),
      rep(minReadCov, nonBlock)
    )
    startPos <- (which(pattern == max(pattern))[1])
    endPos <- which(pattern == max(pattern))[length(which(pattern ==
      max(pattern)))]
  }
  if (fullLeftRight == "Left") {
    pattern <- c(
      rep(cov, blockLength),
      rep(minReadCov, nonBlock)
    )
    startPos <- 1
    endPos <- (which(pattern == min(pattern))[1]) - 1
  }
  if (fullLeftRight == "Right") {
    pattern <- c(
      rep(minReadCov, nonBlock),
      rep(cov, blockLength)
    )
    startPos <- (which(pattern == max(pattern))[1])
    endPos <- length(pattern)
  }
  diff <- mean(abs(viralSubset[, 2] - pattern))
  return(list(
    list(
      diff,
      minReadCov,
      cov,
      "NA",
      startPos,
      endPos,
      "Prophage-like"
    ),
    pattern
  ))
}

#' Full block pattern-translator
#'
#' Translates full block-pattern across a contig. Translate the pattern 1000 bp
#' at a time. Stop translating when the pattern is 5000 bp from the end of the
#' contig.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param bestMatchInfo The information associated with the current best
#'   pattern-match.
#' @param windowSize The window size used to re-average read coverage pileups
#' @param pattern A vector containing the values associated with the block
#'   pattern
#' @return List
#' @keywords internal
blockTranslator <-
  function(viralSubset,
           bestMatchInfo,
           windowSize,
           pattern) {
    minPattern <- min(pattern)
    maxPattern <- max(pattern)
    repeat {
      pattern <- c(
        rep(minPattern, (1000 / windowSize)),
        pattern[-c((length(pattern) - ((
          1000 / windowSize) - 1)):length(pattern))]
      )
      if (pattern[length(pattern) - (5000 / windowSize)] > minPattern) {
        break
      }
      diff <- mean(abs(viralSubset[, 2] - pattern))
      startPos <- which(pattern == max(pattern))[1]
      endPos <- which(pattern == max(pattern))[length(which(pattern ==
        max(pattern)))]
      if (diff < bestMatchInfo[[1]]) {
        bestMatchInfo <- list(
          diff,
          minPattern,
          maxPattern,
          "NA",
          startPos,
          endPos,
          "Prophage-like"
        )
      }
    }
    return(bestMatchInfo)
  }


#' Translate left and right block patterns across contig
#'
#' Translates left and right block patterns across contigs 1000 bp at a time
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param pattern The pattern vector being translated
#' @param leftOrRight Is the left or right block pattern being translated
#' @param windowSize The window size used to re-average read coverage pileups
#' @param minReadCov The baseline value used for the region outside of the block
#'   pattern (either 0 or the minimum VLP-fraction read coverage for the contig)
#' @param cov The height value currently being used for the block pattern
#' @param bestMatchInfo The information associated with the current best
#'   pattern-match.
#' @param minBlockSize The minimum size of the Prophage-like block pattern.
#'   Default is 10,000 bp.
#' @return List
#' @keywords internal
leftRightBlockTranslater <-
  function(viralSubset,
           pattern,
           leftOrRight,
           windowSize,
           minReadCov,
           cov,
           bestMatchInfo,
           minBlockSize) {
    repeat {
      if (leftOrRight == "Left") {
        pattern <- c(
          pattern[-(seq_len(1000 / windowSize))],
          rep(minReadCov, (1000 / windowSize))
        )
        startPos <- 1
        endPos <- (which(pattern == min(pattern))[1]) - 1
      } else {
        pattern <- c(
          rep(minReadCov, (1000 / windowSize)),
          pattern[-c(((length(pattern)) -
            ((
              1000 / windowSize
            ) - 1)):length(pattern))]
        )
        startPos <- (which(pattern == max(pattern))[1])
        endPos <- length(pattern)
      }
      if (length(which(pattern == cov)) <
        (minBlockSize / windowSize) + 1) {
        break
      }
      diff <- mean(abs(viralSubset[, 2] - pattern))
      if (diff < bestMatchInfo[[1]]) {
        bestMatchInfo <- list(
          diff,
          minReadCov,
          cov,
          "NA",
          startPos,
          endPos,
          "Prophage-like"
        )
      }
    }
    return(bestMatchInfo)
  }
