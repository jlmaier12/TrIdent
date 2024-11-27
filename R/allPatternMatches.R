#' Collects pattern-match information for all classifications
#'
#' Collects pattern information associated with all contigs classified as Prophage-like, Sloping and HighCovNoPattern.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @param classifSummTable Classification summary table with whole-community:VLP-fraction read coverage ratios calculated.
#' @keywords internal
allPatternMatches <- function(bestMatchList, classifSummTable){
  A <- 1
  patternMatchList <- list()
  lapply(seq_along(bestMatchList), function(i) {
    if(classifSummTable[which(classifSummTable[,1] == bestMatchList[[i]][[9]]),2] == "NoPattern") return(NULL)
    patternMatchList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(patternMatchList)
}
