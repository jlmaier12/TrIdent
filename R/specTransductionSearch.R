#' Specialized transduction search and plot
#'
#' Search contigs classified as prophage-like for potential specialized
#' transduction and return the plot visualizing the search results.
#'
#' @param VLPpileup A table containing contig names, coverages averaged over 100
#'   bp windows, and contig positions associated with mapping VLP-fraction reads
#'   to whole-community contigs
#' @param classifPatternMatches The pattern match information associated with
#'   each contig classified as prophage-like,  sloping, or HighCovNoPattern
#' @param windowSize The window size used to re-average read coverage pileups
#' @param i The index for the contig currently being assessed
#' @param noReadCov How many bp of no read coverage are encountered before
#'   searching stops? Default is 500.
#' @param specTransLength How many bp of read coverage to look for outside of
#'   prophage borders? Default is 2000.
#' @param classifSumm The summary information associated with each contig
#'   classified as Prophage-like, Sloping, or HighCovNoPattern
#' @param contigName The reference name of the contig currently being assessed
#'   (i.e "NODE_1")
#' @param logScale If TRUE, coverage is plotted in log10. If FALSE, raw coverage
#'   values are plotted. Default is FALSE.
#' @return List containing two objects
#' @keywords internal
specTransductionSearch <- function(contigName,
                                   VLPpileup,
                                   classifPatternMatches,
                                   classifSumm,
                                   windowSize,
                                   i,
                                   noReadCov,
                                   specTransLength,
                                   logScale) {
    specTransSumm <- c(contigName, rep(NA, 5))
    viralSubsetZoom <- prophageLikeZoom(
        VLPpileup[which(VLPpileup[, 1] ==
            contigName), ],
        classifPatternMatches, i, 500, windowSize
    )
    borders <-
        prophageLikeBorders(
            VLPpileup[which(VLPpileup[, 1] == contigName), ],
            classifPatternMatches, i, windowSize
        )
    startPosBp <- borders[[1]]
    startPosRow <- which(viralSubsetZoom[, 3] == startPosBp)
    endPosBp <- borders[[2]]
    endPosRow <- which(viralSubsetZoom[, 3] == endPosBp)
    viralSubsetZoom$logcoverage <- abs(log10(viralSubsetZoom[, 2]))
    viralSubsetZoom[viralSubsetZoom == Inf] <- 0
    startSearch <-
        zeroCountSearch("start", viralSubsetZoom, startPosRow, noReadCov)
    endSearch <-
        zeroCountSearch("end", viralSubsetZoom, endPosRow, noReadCov)
    if (startSearch[[2]] - startSearch[[1]] >= specTransLength / 100) {
        SpecTransLeft <- viralSubsetZoom[startPosRow - (startSearch[[2]] -
            startSearch[[1]]), 3]
        specTransSumm[c(3, 5)] <- c("yes", (startSearch[[2]] -
            startSearch[[1]]) * 100)
    } else {
        SpecTransLeft <- startPosBp
        specTransSumm[c(3, 5)] <- c("no", NA)
    }
    if (endSearch[[2]] - endSearch[[1]] >= specTransLength / 100) {
        specTransRight <- viralSubsetZoom[endPosRow + (endSearch[[2]] -
            endSearch[[1]]), 3]
        specTransSumm[c(4, 6)] <-
            c("yes", (endSearch[[2]] - endSearch[[1]]) * 100)
    } else {
        specTransRight <- endPosBp
        specTransSumm[c(4, 6)] <- c("no", NA)
    }
    specTransSumm[2] <-
        ifelse((
            startSearch[[2]] - startSearch[[1]] >= specTransLength / 100 |
                endSearch[[2]] - endSearch[[1]] >= specTransLength / 100
        ),
        "yes",
        "no"
        )
    plot <- specTransductionPlot(
        viralSubsetZoom,
        startPosBp,
        endPosBp,
        SpecTransLeft,
        specTransRight,
        contigName,
        classifPatternMatches,
        i,
        specTransSumm,
        logScale,
        classifSumm
    )
    return(list(specTransSumm, plot))
}
