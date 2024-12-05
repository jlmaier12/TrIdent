#' Change the read coverage rolling mean window size
#'
#' Re-averages window sizes of read coverage averages. Start with 100bp windows always. Cannot make window size less than 100bp.
#'
#' @param cleanPileup A read coverage dataset that has been cleaned and reformatted by the readcovdf_formatter function
#' @param windowSize The number of base pairs to average coverage values over. Options are 100, 500, 1000, or 2000 only!
#' @keywords internal
changeWindowSize <- function(cleanPileup, windowSize){
  coverage <- rep(NA, nrow(cleanPileup))
  X <- 0
  Y <- windowSize / 100
  A <- 1
  repeat{
    coverage[A] <- mean(cleanPileup[c(X:Y),2])
    A <- A + 1
    X <- X + (windowSize / 100)
    Y <- Y + (windowSize / 100)
    if (Y > nrow(cleanPileup)) break
  }
  coverage <- coverage[!is.na(coverage)]
  position <- seq(windowSize, length(coverage) * windowSize, windowSize)
  contigName <- rep(cleanPileup[1,1], length(position))
  newWindowSizePileup <- cbind.data.frame(contigName, coverage, position) %>% as.data.frame()
  newWindowSizePileup[NARemover(newWindowSizePileup)] <- 0
  return(newWindowSizePileup)
}
