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
  lys <- lapply(seq_along(bestMatchList), function(i) {
    if (classifSummTable[which(classifSummTable[, 1] ==
                               bestMatchList[[i]][[8]]), 2] == "NoPattern") {
      return(NULL)
    }
    bestMatchList[[i]]
  })
  return(lys[!vapply(lys, is.null, logical(1))])
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
  lys <- lapply(seq_along(bestMatchList), function(i) {
    if (bestMatchList[[i]][[7]] != "Prophage-like") {
      return(NULL)
    }
    bestMatchList[[i]]
  })
  return(lys[!vapply(lys, is.null, logical(1))])
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
  lys <- lapply(seq_along(bestMatchList), function(i) {
    if (bestMatchList[[i]][[7]] != "Sloping") {
      return(NULL)
    }
    bestMatchList[[i]]
  })
  return(lys[!vapply(lys, is.null, logical(1))])
}
