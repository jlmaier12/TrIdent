#' Change slope of sloping pattern
#'
#' Change the value of the slope used for the sloping pattern-match
#'
#' @param leftOrRight
#'  Generate pattern for negative slope
#'  (left to right, i.e. 'Left') or
#'  positive slope (right to left, i.e. 'Right')
#' @param slopeBottom
#'  The value for the bottom of the sloping value
#' @param halfToMaxReadCov
#'  Half of the max VLP-fraction read coverage
#'  divided by 10
#' @param cov
#'  The value for the top of the slope
#' @param viralSubset
#'  A subset of the read coverage pileup that pertains
#'  only to the contig currently being assessed
#' @param windowSize
#'  The window size used to re-average read coverage pileup
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
        slope <- ifelse(
            leftOrRight == "Left",
            (slopeBottom - cov) / (nrow(viralSubset) * windowSize),
            (cov - slopeBottom) / (nrow(viralSubset) * windowSize)
        )
        return(list(
            diff,
            slopeBottom,
            cov,
            covSteps,
            1,
            length(pattern),
            slope,
            "Sloping"
        ))
    }
