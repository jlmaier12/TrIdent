#' No pattern pattern-match
#'
#' A horizontal line at the mean or median coverage should be an optimal
#' pattern-match if the contig read coverage displays no sloping or block
#' patterns
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @return List
#' @keywords internal
noPattern <- function(viralSubset) {
    pattern1 <- rep(median(viralSubset[, 2]), nrow(viralSubset))
    pattern2 <- rep(mean(viralSubset[, 2]), nrow(viralSubset))
    diff1 <- mean(abs(viralSubset[, 2] - pattern1))
    diff2 <- mean(abs(viralSubset[, 2] - pattern2))
    diff <- ifelse(diff1 < diff2, diff1, diff2)
    value <- ifelse(diff1 < diff2, median(viralSubset[, 2]), mean(viralSubset[, 2]))
    bestMatchInfo <-
        list(
            diff,
            value,
            nrow(viralSubset),
            "NA",
            1,
            nrow(viralSubset),
            "NoPattern"
        )
    return(bestMatchInfo)
}
