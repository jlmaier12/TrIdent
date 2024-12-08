#' Collects Prophage-like classification pattern-match information
#'
#' Collects pattern information associated with all contigs classified as Prophage-like.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @return List
#' @keywords internal
allProphageLikeClassifs <- function(bestMatchList){
  A <- 1
  prophageLikeClassifList <- list()
  lapply(seq_along(bestMatchList), function(i) {
    classification <-  bestMatchList[[i]][[8]]
    if(classification != "Prophage-like") return(NULL)
    prophageLikeClassifList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(prophageLikeClassifList)
}
