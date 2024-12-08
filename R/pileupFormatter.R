#' Correctly formats pileup files.
#'
#' Places columns in correct order and renames columns. Cleans the contig labels to remove excess information after whitespace.
#'
#' @param pileup A table containing contig names, read coverages averaged over 100 bp windows, and contig positions
#' @return dataframe
#' @keywords internal
pileupFormatter <- function(pileup) {
  colClasses <- lapply(seq_along(pileup), function(i) class(pileup[,i]))
  for (i in c(which(colClasses == "integer"))){
    if (length(which(pileup[,i] == 100)) > 1) posColIdx <- i
  }
  cleanPileup <- cbind.data.frame(pileup[,which(colClasses == "character")], pileup[,which(colClasses == "numeric")], pileup[,posColIdx])
  colnames(cleanPileup) <- c("contigName", "coverage", "position")
  cleanPileup$contigName <- gsub("\\s.*", "", cleanPileup$contigName)
  return(cleanPileup)
}
