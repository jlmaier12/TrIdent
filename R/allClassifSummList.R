#' Collects Prophage-like and Sloping classification pattern-match information
#'
#' Collects pattern-match information associated with all contigs classified as Prophage-like or Sloping.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @keywords internal
allClassifSummList <- function(bestMatchList){
  A <- 1
  classifList <- list()
  if (length(bestMatchList) == 0) stop ("NO TRANSDUCTION EVENTS FOUND")
  lapply(seq_along(bestMatchList), function(i) {
    classification <- bestMatchList[[i]][[8]]
    if(classification == "NoPattern") return(NULL)
    classifList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(classifList)
}
