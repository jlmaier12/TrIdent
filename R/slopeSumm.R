#' Summarize slopes for sloping classifications
#'
#' Add slope information for sloping classifications to summary table
#'
#' @param classifSumm Classification summary table
#' @param slopingClassifList A list containing pattern match information
#'   associated with all contigs classified as sloping.
#' @param windowSize The window size used to re-average read coverage pileups
#' @return dataframe
#' @keywords internal
slopeSumm <- function(classifSumm, slopingClassifList, windowSize) {
    classifSumm$slope <- rep(NA, nrow(classifSumm))
    if (length(slopingClassifList) == 0) {
        return(classifSumm)
    }
    lapply(seq_along(slopingClassifList), function(i) {
        classifSumm[which(classifSumm[, 1] ==
            slopingClassifList[[i]][[8]]), 10] <<-
            round(slopingClassifList[[i]][[4]] / windowSize, digits = 4)
    })
    return(classifSumm)
}
