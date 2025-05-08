#' Create histogram of normalized pattern-match scores
#'
#' Plots a histogram of normalized match scores for all Prophage-like, Sloping
#' and HighCovNoPattern classifications and colors the plot based on the
#' classifications. 
#'
#' @param summaryList Classification summary table filtered to only include
#'   contigs with Prophage-like, Sloping and HighCovNoPattern classifications
#' @return ggplot object
#' @keywords internal
resultsHisto <- function(summaryList) {
  normMatchScore <- classifications <- NULL
  summaryTable <- summaryList[[1]]
  cleanSummaryTable <- summaryList[[2]]
  brks <- seq(
    min(cleanSummaryTable$normMatchScore),
    max(cleanSummaryTable$normMatchScore),
    length.out = 40
  )
    plot <- ggplot(data = cleanSummaryTable) +
      theme_bw() +
      geom_histogram(aes(x = normMatchScore, fill = classifications),
        breaks = brks
      ) +
      labs(
        title = "Distribution of Pattern-Match Classifications",
        x = "Normalized Pattern-Match Score",
        y = "count",
        caption = paste(
          "(Lower scores are better pattern-matches) \n",
          (length(which(summaryTable[, 2] == "NoPattern"))),
          "contigs were classified as 'NoPattern' and are not
                    displayed on this plot."
        )
      ) +
      theme(
        plot.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.caption = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
  return(list(plot))
}
