#' Determine Prophage-like read coverage elevation in whole-community
#'
#' Determines whether a detected Prophage-like genetic element has read coverage in the whole-community
#' that is either elevated or depressed compared to the average read coverage of the non-prophage region.
#'
#' @param classifSummTable Classification summary table
#' @param prophageLikeClassifList A list containing pattern match information associated with all contigs classified as Prophage-like.
#' @param VLPpileup A table containing contig names, coverages averaged over 100 bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @param WCpileup A table containing contig names, coverages averaged over 100 bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param windowSize The window size used to re-average read coverage pileups
#' @keywords internal
prophageLikeElevation <- function(classifSummTable, prophageLikeClassifList, VLPpileup, WCpileup, windowSize){
  message("Identifying highly active/abundant or heterogenously integrated prophage-like elements")
  classifSummTable$prophageLikeRegionReadCov <- rep(NA, nrow(classifSummTable))
  classifSummTable$prophageLikeElevationRatio <- rep(NA, nrow(classifSummTable))
  if(length(prophageLikeClassifList) == 0) return(classifSummTable)
  lapply(seq_along(prophageLikeClassifList), function(i) {
    viralSubset <- changeWindowSize(VLPpileup[which(VLPpileup[,1] == prophageLikeClassifList[[i]][[9]]),], windowSize)
    startPos <- prophageLikeClassifList[[i]][[5]]
    endPos <- prophageLikeClassifList[[i]][[6]]
    contigName <- prophageLikeClassifList[[i]][[9]]
    blockLengthBp <- abs(endPos - startPos) * windowSize
    nonBlockLengthBp <- (nrow(viralSubset) * windowSize) - blockLengthBp
    if (nonBlockLengthBp < 20000) classifSummTable[which(classifSummTable[,1] == contigName), c(7:8)] <<- "CBD" ##can't be determined
    else {
      microbialSubset <- changeWindowSize(WCpileup[which(WCpileup[, 1] == prophageLikeClassifList[[i]][[9]]),], windowSize)
      prophageLikeRegion <- microbialSubset[c(startPos:endPos), 2]
      nonProphageLikeRegion <- microbialSubset[which(!microbialSubset[,2] %in% prophageLikeRegion),2]
      ratio <- round(mean(prophageLikeRegion) / mean(nonProphageLikeRegion), digits=4)
      classifSummTable[which(classifSummTable[,1] == contigName), 8] <<- ratio
      if(ratio > 1.3) classifSummTable[which(classifSummTable[,1] == contigName),7] <<- "Elevated"
      else if (ratio < 0.75) classifSummTable[which(classifSummTable[,1] == contigName),7] <<- "Depressed"
      else classifSummTable[which(classifSummTable[,1] == contigName),7] <<- "None"
    }
  })
  return(classifSummTable)
}
