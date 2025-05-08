#' Change the read coverage rolling mean window size
#'
#' Re-averages window sizes of read coverage averages. Start with 100bp windows
#' always. Cannot make window size less than 100bp.
#'
#' @param cleanPileup A read coverage dataset that has been cleaned and
#'   reformatted.
#' @param windowSize The number of base pairs to average coverage values over.
#'   Options are 100, 500, 1000, or 2000 only!
#' @return Dataframe
#' @keywords internal
changeWindowSize <- function(cleanPileup, windowSize) {
  coverage <- rep(NA, nrow(cleanPileup))
  X <- 0
  Y <- windowSize / 100
  A <- 1
  repeat {
    coverage[A] <- mean(cleanPileup[c(X:Y), 2])
    A <- A + 1
    X <- X + (windowSize / 100)
    Y <- Y + (windowSize / 100)
    if (Y > nrow(cleanPileup)) {
      break
    }
  }
  coverage <- coverage[!is.na(coverage)]
  position <-
    seq(windowSize, length(coverage) * windowSize, windowSize)
  contigName <- rep(cleanPileup[1, 1], length(position))
  newWindowSizePileup <- cbind.data.frame(
    contigName, coverage, position)
  newWindowSizePileup[NARemover(newWindowSizePileup)] <- 0
  return(newWindowSizePileup)
}

#' Correctly formats pileup files.
#'
#' Places columns in correct order and renames columns. Cleans the contig labels
#' to remove excess information after whitespace.
#'
#' @param pileup A table containing contig names, read coverages averaged over
#'   100 bp windows,and contig positions
#' @return dataframe
#' @keywords internal
pileupFormatter <- function(pileup) {
  colClasses <- vapply(pileup, class, character(1))
  for (i in c(which(colClasses == "integer"))) {
    if (length(which(pileup[, i] == 100)) > 1) {
      posColIdx <- i
    }
  }
  cleanPileup <-
    cbind.data.frame(
      pileup[, which(colClasses == "character")],
      pileup[, which(colClasses == "numeric")],
      pileup[, posColIdx]
    )
  colnames(cleanPileup) <- c("contigName", "coverage", "position")
  cleanPileup$contigName <- gsub("\\s.*", "", cleanPileup$contigName)
  return(cleanPileup)
}

#' NA remover
#'
#' Removes NAs from dataframe.
#'
#' @seealso
#' \url{https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame/18143097#18143097}
#'
#' @param x dataset with potential NAs
#' @return Dataframe
#' @keywords internal
NARemover <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}
