#' Whole-community:VLP-fraction read coverage ratio calculator
#'
#' Calculate the whole-community:VLP-fraction read coverage ratio for every contig using the median
#' read coverage values. If the ratio is less than 2 (i.e VLP-fraction read coverage
#' is, on average, at least half the whole-community read coverage), then the contig is classified as HighCovNoPattern
#'
#' @param classifSummTable Classification summary table
#' @param WCpileup A table containing contig names, coverages averaged over 100 bp windows, and contig positions associated with mapping whole-community reads to whole-community contigs
#' @param VLPpileup  A table containing contig names, coverages averaged over 100 bp windows, and contig positions associated with mapping VLP-fraction reads to whole-community contigs
#' @keywords internal
microbialToViralRatioCalc <- function(classifSummTable, WCpileup, VLPpileup){
  noneClassifIdxs <- which(classifSummTable[,2] == "NoPattern")
  if(length(noneClassifIdxs) == 0) return(classifSummTable)
  lapply(seq_along(noneClassifIdxs), function(p) {
    i <- noneClassifIdxs[[p]]
    contigName <- classifSummTable[i,1]
    viralSubset <- VLPpileup[which(VLPpileup[,1] == contigName),]
    viralSubset[NARemover(viralSubset)] <- 0
    microbialSubset <- WCpileup[which(WCpileup[,1] == contigName),]
    microbialSubset[NARemover(microbialSubset)] <- 0
    microbialToViralRatio <- median(microbialSubset[,2]) / median(viralSubset[,2])
    classifSummTable[i,2] <<- ifelse(microbialToViralRatio < 2, "HighCovNoPattern", "NoPattern")
  })
  return(classifSummTable)
}
