#' Collects pattern-match information for all classifications
#'
#' Collects pattern information associated with all contigs classified as
#' Prophage-like, Sloping and HighCovNoPattern.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @param classifSummTable Classification summary table with
#'   whole-community:VLP-fraction read coverage ratios calculated.
#' @return List
#' @keywords internal
allPatternMatches <- function(bestMatchList, classifSummTable) {
  A <- 1
  patternMatchList <- list()
  lapply(seq_along(bestMatchList), function(i) {
    if (classifSummTable[which(classifSummTable[, 1] ==
      bestMatchList[[i]][[8]]), 2] == "NoPattern") {
      return(NULL)
    }
    patternMatchList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(patternMatchList)
}


#' Collects Prophage-like classification pattern-match information
#'
#' Collects pattern information associated with all contigs classified as
#' Prophage-like.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @return List
#' @keywords internal
allProphageLikeClassifs <- function(bestMatchList) {
  A <- 1
  prophageLikeClassifList <- list()
  lapply(seq_along(bestMatchList), function(i) {
    classification <- bestMatchList[[i]][[7]]
    if (classification != "Prophage-like") {
      return(NULL)
    }
    prophageLikeClassifList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(prophageLikeClassifList)
}

#' Collects Sloping classification pattern-match i nformation
#'
#' Collects pattern information associated with all contigs classified as
#' Sloping in the patternMatcher function.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @return List
#' @keywords internal
allSlopingClassifs <- function(bestMatchList) {
  A <- 1
  slopingClassifList <- list()
  lapply(seq_along(bestMatchList), function(i) {
    classification <- bestMatchList[[i]][[7]]
    if (classification != "Sloping") {
      return(NULL)
    }
    slopingClassifList[[A]] <<- bestMatchList[[i]]
    A <<- A + 1
  })
  return(slopingClassifList)
}
