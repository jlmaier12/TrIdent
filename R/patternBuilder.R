#' Pattern-builder
#'
#' Builds the pattern (vector) associated with the best pattern-match' for each
#' contig classified as Prophage-like, Sloping, or HighCovNoPattern.
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param classifList A list containing pattern match information associated
#'   with all classified contigs.
#' @param i The list index associated with each contig's pattern-match
#'   information
#' @param classification The contig's classification assigned by the
#'   TrIdentClassifier function
#' @return Vector
#' @keywords internal
patternBuilder <-
    function(viralSubset,
             classifList,
             classification,
             i) {
        minReadCov <- classifList[[i]][[2]]
        maxReadCov <- classifList[[i]][[3]]
        slopingCovSteps <- classifList[[i]][[4]]
        startPos <- classifList[[i]][[5]]
        endPos <- classifList[[i]][[6]]
        if (classification == "Prophage-like") {
            if (startPos == 1) {
                pattern <- c(
                    rep(maxReadCov, endPos),
                    rep(
                        minReadCov,
                        (nrow(viralSubset) - endPos)
                    )
                )
            } else if (endPos == nrow(viralSubset)) {
                pattern <- c(
                    rep(
                        minReadCov,
                        startPos
                    ),
                    rep(
                        maxReadCov,
                        (nrow(viralSubset) -
                            startPos)
                    )
                )
            } else {
                blockRegion <- endPos - startPos
                pattern <-
                    c(
                        rep(minReadCov, startPos),
                        rep(maxReadCov, blockRegion),
                        rep(
                            minReadCov,
                            (
                                nrow(viralSubset) - (blockRegion + startPos)
                            )
                        )
                    )
            }
        } else if (classification == "Sloping") {
            if (startPos == 1 &
                endPos == nrow(viralSubset) &
                slopingCovSteps < 0) {
                pattern <- seq(maxReadCov, minReadCov, slopingCovSteps)
            } else if (startPos == 1 &
                endPos == nrow(viralSubset) &
                slopingCovSteps > 0) {
                pattern <- seq(minReadCov, maxReadCov, slopingCovSteps)
            } else if (startPos != 1) {
                pattern <- c(
                    rep(
                        min(viralSubset[, 2]),
                        startPos - 1
                    ),
                    seq(
                        maxReadCov,
                        minReadCov,
                        slopingCovSteps
                    )
                )
            } else if (endPos != nrow(viralSubset)) {
                pattern <- c(
                    seq(
                        minReadCov,
                        maxReadCov,
                        slopingCovSteps
                    ),
                    rep(
                        min(viralSubset[, 2]),
                        (nrow(viralSubset) -
                            endPos)
                    )
                )
            }
            ## for NoPattern, min read cov = med read cov and
            ## max read cov = contig length
        } else if (classification == "HighCovNoPattern") {
            pattern <- rep(
                minReadCov,
                maxReadCov
            )
        }
        return(pattern)
    }
