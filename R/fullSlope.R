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
                if (abs(slopeChangeInfoL[[7]]) <
                    minSlope | slopeChangeInfoL[[7]] > 0) {
                    break
                }
                if (abs(slopeChangeInfoR[[7]]) <
                    minSlope | slopeChangeInfoR[[7]] < 0) {
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
