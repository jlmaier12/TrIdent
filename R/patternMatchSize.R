#' Pattern-match size calculator
#'
#' Calculate the size (bp) of the matching region for Prophage-like and
#' Sloping patterns
#'
#' @param classifSumm Classification summary table
#' @param classifList
#' A list containing pattern match information associated with all contig
#' classifications
#' @param windowSize The window size used to re-average read coverage pileups
#' @return dataframe
#' @keywords internal
patternMatchSize <- function(classifSumm, classifList, windowSize){
message("Determining sizes (bp) of pattern matches")
classifSumm <- as.data.frame(classifSumm)
classifSumm$matchSize <- rep(NA, nrow(classifSumm))
classifSumm$startPosBp <- rep(NA, nrow(classifSumm))
classifSumm$endPosBp <- rep(NA, nrow(classifSumm))
lapply(seq_along(classifList), function(i) {
    contigName <- classifList[[i]][[9]]
    startPos <- classifList[[i]][[5]]
    endPos <- classifList[[i]][[6]]
    classifSumm[which(classifSumm[,1] == contigName),5] <<- (length(c(startPos:endPos)) - 1) * windowSize
    classifSumm[which(classifSumm[,1] == contigName),6] <<- startPos * windowSize
    classifSumm[which(classifSumm[,1] == contigName),7] <<- endPos * windowSize
  })
return(classifSumm)
}
