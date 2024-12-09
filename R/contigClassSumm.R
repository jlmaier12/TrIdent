#' Summarizes pattern-match information
#'
#' Summarizes the classifications made in the patternMatcher() function into a
#' dataframe.
#'
#' @param bestMatchList Classifications made with patternMatcher function.
#' @return dataframe
#' @keywords internal
contigClassSumm <- function(bestMatchList) {
    if (length(bestMatchList) == 0) {
        stop("NO TRANSDUCTION EVENTS FOUND")
    }
    contigName <- rep(NA, length(bestMatchList))
    classifications <- rep(NA, length(bestMatchList))
    normMatchScore <- rep(NA, length(bestMatchList))
    lapply(seq_along(bestMatchList), function(i) {
        normMatchScore[i] <<- bestMatchList[[i]][[10]]
        contigName[i] <<- bestMatchList[[i]][[9]]
        classifications[i] <<- bestMatchList[[i]][[8]]
    })
    classifSumm <-
        cbind.data.frame(contigName, classifications, normMatchScore)
    return(classifSumm)
}
