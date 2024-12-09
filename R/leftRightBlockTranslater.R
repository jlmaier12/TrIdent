#' Translate left and right block patterns across contig
#'
#' Translates left and right block patterns across contigs 1000 bp at a time
#'
#' @param viralSubset
#' A subset of the read coverage pileup that pertains only to the contig
#' currently being assessed
#' @param pattern The pattern vector being translated
#' @param leftOrRight Is the left or right block pattern being translated
#' @param windowSize The window size used to re-average read coverage pileups
#' @param minReadCov
#' The baseline value used for the region outside of the block pattern (either
#' 0 or the minimum VLP-fraction read coverage for the contig)
#' @param cov The height value currently being used for the block pattern
#' @param bestMatchInfo
#' The information associated with the current best pattern match.
#' @param minBlockSize
#' The minimum size of the Prophage-like block pattern. Default is 10,000 bp.
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
            pattern <- if (leftOrRight == "Left") {
                c(
                    pattern[-(seq_len(1000 / windowSize))],
                    rep(minReadCov, (1000 / windowSize))
                )
            } else {
                c(
                    rep(minReadCov, (1000 / windowSize)),
                    pattern[-c(((length(pattern)) -
                        ((
                            1000 / windowSize
                        ) - 1)):length(pattern))]
                )
            }
            if (length(which(pattern == cov)) < (minBlockSize / windowSize) + 1) {
                break
            }
            diff <- mean(abs(viralSubset[, 2] - pattern))
            startPos <- ifelse(leftOrRight == "Left",
                1, (which(pattern == max(pattern))[1])
            )
            endPos <- ifelse(leftOrRight == "Left",
                (which(pattern == min(pattern))[1]) - 1, length(pattern)
            )
            if (diff < bestMatchInfo[[1]]) {
                bestMatchInfo <- list(
                    diff,
                    minReadCov,
                    cov,
                    "NA",
                    startPos,
                    endPos,
                    "NA",
                    "Prophage-like"
                )
            }
        }
        return(bestMatchInfo)
    }
