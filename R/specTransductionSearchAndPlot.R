#' Specialized transduction search and plot
#'
#' Search contigs classified as prophage-like for potential specialized
#' transduction and return the plot visualizing the search results.
#'
#' @param VLPpileup
#' A table containing contig names, coverages averaged over 100 bp windows, and
#' contig positions associated with mapping VLP-fraction reads to
#' whole-community contigs
#' @param classifPatternMatches
#' The pattern match information associated with each contig classified as
#' prophage-like,  sloping, or HighCovNoPattern
#' @param windowSize The window size used to re-average read coverage pileups
#' @param i The index for the contig currently being assessed
#' @param noReadCov
#' How many bp of no read coverage are encountered before searching stops?
#' Default is 500.
#' @param specTransLength
#' How many bp of read coverage to look for outside of prophage borders?
#' Default is 2000.
#' @param classifSumm
#' The summary information associated with each contig classified as
#' Prophage-like, Sloping, or HighCovNoPattern
#' @param contigName
#' The reference name of the contig currently being assessed (i.e "NODE_1")
#' @param logScale
#' If TRUE, coverage is plotted in log10. If FALSE, raw coverage values are
#' plotted. Default is FALSE.
#' @return List containing two objects
#' @keywords internal
specTransductionSearchAndPlot <- function(contigName, VLPpileup,
                                          classifPatternMatches, classifSumm,
                                          windowSize, i, noReadCov,
                                          specTransLength, logScale){
position <- logcoverage <- NULL
specTransSumm <- c(contigName, rep(NA, 5))
if (classifSumm[which(classifSumm[,1] == contigName),6] == "Elevated") prophageLikeInfo <- "Highly active/abundant prophage-like element"
else if (classifSumm[which(classifSumm[,1] == contigName),6] == "Depressed") prophageLikeInfo <- "Not homogenously present/integrated prophage-like element"
else prophageLikeInfo <- NULL
viralSubsetZoom <- prophageLikeZoom(VLPpileup[which(VLPpileup[,1] == contigName),],
                                    classifPatternMatches, i, 500, windowSize)
borders <- prophageLikeBorders(VLPpileup[which(VLPpileup[,1] == contigName),],
                               classifPatternMatches, i, windowSize)
startPosBp <- borders[[1]]
startPosRow <- which(viralSubsetZoom[,3] == startPosBp)
endPosBp <- borders[[2]]
endPosRow <- which(viralSubsetZoom[,3] == endPosBp)
viralSubsetZoom$logcoverage <- abs(log10(viralSubsetZoom[,2]))
viralSubsetZoom[viralSubsetZoom == Inf] <- 0
startSearch <- zeroCountSearch("start", viralSubsetZoom, startPosRow, noReadCov)
endSearch <- zeroCountSearch("end", viralSubsetZoom, endPosRow, noReadCov)
if (startSearch[[2]] - startSearch[[1]] >= specTransLength / 100) {
    SpecTransLeft <- viralSubsetZoom[startPosRow - (startSearch[[2]] - startSearch[[1]]),3]
    specTransSumm[c(3,5)] <- c("yes", (startSearch[[2]] - startSearch[[1]]) * 100)
} else {
    SpecTransLeft <- startPosBp
    specTransSumm[c(3,5)] <- c("no", NA)
}
if (endSearch[[2]] - endSearch[[1]] >= specTransLength / 100) {
    specTransRight <- viralSubsetZoom[endPosRow + (endSearch[[2]] - endSearch[[1]]),3]
    specTransSumm[c(4,6)] <- c("yes", (endSearch[[2]] - endSearch[[1]]) * 100)
} else {
    specTransRight <- endPosBp
    specTransSumm[c(4,6)] <- c("no", NA)
  }
specTransSumm[2] <- ifelse((startSearch[[2]] - startSearch[[1]] >= specTransLength / 100 | endSearch[[2]] - endSearch[[1]] >= specTransLength / 100),
                             "yes", "no")
fill <- ifelse(specTransSumm[2] == "yes", "seagreen", "deepskyblue3")
alphaL <- ifelse(specTransSumm[3] == "yes", 1, 0)
alphaR <- ifelse(specTransSumm[4] == "yes", 1, 0)
coverageType <- if(logScale == TRUE) viralSubsetZoom$logcoverage else viralSubsetZoom$coverage
plot <- (ggplot(data=as.data.frame(viralSubsetZoom), aes(x=position, y=coverageType))+
             geom_area(fill=fill) +
             geom_vline(xintercept=c(startPosBp, endPosBp), linewidth=1)+
             geom_vline(xintercept=SpecTransLeft, color="red", alpha=alphaL,
                        linewidth=1)+
             geom_vline(xintercept=specTransRight, color="red", alpha=alphaR,
                        linewidth=1)+
             scale_x_continuous(expand = c(0, 0)) +
             theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   plot.subtitle=element_text(size=10),
                   panel.background=element_blank(),
                   axis.line=element_line(colour="black"),
                   text=element_text(size=15))+
             labs(title=paste(contigName, classifPatternMatches[[i]][[8]],
                              prophageLikeInfo),
                  subtitle=paste("Specialized transduction on left:",
                                 specTransSumm[3],
                                 ",", "on right:", specTransSumm[4]),
                  x="Contig Position (bp)",
                  y=paste("VLP-Fraction Read Coverage", ifelse(logScale == TRUE,
                                                               "(Log10)", ""))))
return(list(specTransSumm, plot))
}
