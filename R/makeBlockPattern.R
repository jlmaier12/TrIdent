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
