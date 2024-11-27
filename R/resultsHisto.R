#' Create histogram of normalized pattern-match scores
#'
#' Plots a histogram of normalized match scores for all Prophage-like, Sloping and HighCovNoPattern classifications
#' and colors the plot based on the classifications. A suggested filtering threshold is provided for filtering results
#' based on the quality of the pattern-match.
#'
#' @param summaryList Classification summary table filtered to only include contigs with Prophage-like, Sloping and HighCovNoPattern classifications
#' @keywords internal
resultsHisto <- function (summaryList){
  normMatchScore <- classifications <- NULL
  summaryTable <- summaryList[[1]]
  cleanSummaryTable <- summaryList[[2]]
  brks <- seq(min(cleanSummaryTable$normMatchScore), max(cleanSummaryTable$normMatchScore), length.out=40)
  histogram <- hist(cleanSummaryTable$normMatchScore, breaks=brks, plot=FALSE)
  ST <- suggThreshold(cleanSummaryTable$normMatchScore, histogram)
  plot <- ggplot(data=cleanSummaryTable)+
    theme_bw()+
    geom_histogram(aes(x=normMatchScore, fill=classifications), breaks=brks)+
    geom_vline(xintercept=ST, color="red")+
    annotate(geom="label", x=ST, y=max(histogram$counts),label=round(ST, digits=2)) +
    labs(title="Quality of pattern matches", x="Normalized Pattern Match Score",
         y="count", caption=paste("(Lower scores are better matches) \n
         Suggested Filtering Threshold=",round(ST, digits=2), "\n",
         (length(which(summaryTable[,2] == "NoPattern"))),
         "contigs were classified as 'NoPattern' and are not displayed on this plot."))
  return(plot)
}
