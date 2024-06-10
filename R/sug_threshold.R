#' Calculates a suggested filtering threshold for plotting/saving TrIdent results
#'
#' Uses the distribution of normalized pattern match-scores to calculate a suggested filtering threshold.
#'
#' @param norm_match_scores Mormalized pattern-match score vector
#' @param histogram A histogram object created from normaized pattern-match scores
#' @keywords internal
sug_threshold <- function (norm_match_scores, histogram){
  max <- which.max(histogram$counts)
  SD<- sd(histogram$counts)/2
  zeros <- c()
  if(max(norm_match_scores) <= 0.3) return(max(norm_match_scores))
  repeat{
    freq1 <- histogram$counts[[max]]
    max <- max+1
    freq2 <- histogram$counts[[max]]
    if(freq1-freq2 == 0) {
      zeros <- c(zeros, max)
    }
    if(freq1-freq2 < 0) break
    if(max==length( histogram$counts)) {
      max <- zeros[[(length(zeros))]]
      break
    }
  }
  QCvalue <- histogram$breaks[ 0:1 + (max+SD) ]
  return(QCvalue[[1]])
}


