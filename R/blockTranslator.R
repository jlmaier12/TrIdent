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
