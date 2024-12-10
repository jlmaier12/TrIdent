#' VLP-fraction:whole-community read coverage ratio calculator
#'
#' Calculate the VLP-fraction:whole-community read coverage ratio for every
#' contig using the median read coverage values. If the ratio is greater than 2
#' (i.e VLP-fraction read coverage is, on average, at least double the
#' whole-community read coverage), then the contig is classified as
#' HighCovNoPattern
#'
#' @param classifSumm Classification summary table
#' @param WCpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping whole-community
#'   reads to whole-community contigs
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping VLP-fraction reads
#'   to whole-community contigs
#' @return dataframe
#' @keywords internal
VLPtoWCRatioCalc <- function(classifSumm, WCpileup, VLPpileup) {
    classifSumm$VLPWCRatio <- rep(NA, nrow(classifSumm))
    noneClassifIdxs <- which(classifSumm[, 2] == "NoPattern")
    if (length(noneClassifIdxs) == 0) {
        return(classifSumm)
    }
    lapply(seq_along(noneClassifIdxs), function(p) {
        i <- noneClassifIdxs[[p]]
        contigName <- classifSumm[i, 1]
        viralSubset <- VLPpileup[which(VLPpileup[, 1] == contigName), ]
        viralSubset[NARemover(viralSubset)] <- 0
        microbialSubset <-
            WCpileup[which(WCpileup[, 1] == contigName), ]
        microbialSubset[NARemover(microbialSubset)] <- 0
        VLPtoWCratio <-
            round(median(viralSubset[, 2] / median(microbialSubset[, 2])),
                digits = 4
            )
        classifSumm[i, 2] <<-
            ifelse(VLPtoWCratio > 2, "HighCovNoPattern", "NoPattern")
        classifSumm[i, 4] <<- VLPtoWCratio
    })
    return(classifSumm)
}
