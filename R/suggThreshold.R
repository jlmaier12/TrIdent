#' Calculates a suggested filtering threshold for plotting/saving TrIdent results
#'
#' Uses the distribution of normalized pattern match-scores to calculate a suggested filtering threshold. The
#' threshold is two standard deviations to the right of the maximum histogram value.
#'
#' @param normMatchScores Normalized pattern-match score vector
#' @param histogram A histogram object created from normalized pattern-match scores
#' @keywords internal
suggThreshold <- function (normMatchScores, histogram){
  max <- which.max(histogram$counts)
  SD <- sd(histogram$counts) * 2
  QCvalue <- histogram$breaks[0:1 + (max + SD)]
  return(QCvalue[[1]])
}


