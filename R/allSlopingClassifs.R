#' Collects Sloping classification pattern-match information
#'
#' Collects pattern information associated with all
#' contigs classified as Sloping in the patternMatcher function.
#'
#' @param bestMatchList
#'  Classifications made with patternMatcher function.
#' @return List
#' @keywords internal
allSlopingClassifs <- function(bestMatchList) {
    A <- 1
    slopingClassifList <- list()
    lapply(seq_along(bestMatchList), function(i) {
        classification <- bestMatchList[[i]][[8]]
        if (classification != "Sloping") {
            return(NULL)
        }
        slopingClassifList[[A]] <<- bestMatchList[[i]]
        A <<- A + 1
    })
    return(slopingClassifList)
}
