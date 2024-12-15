#' Create histogram of normalized pattern-match scores
#'
#' Plots a histogram of normalized match scores for all Prophage-like, Sloping
#' and HighCovNoPattern classifications and colors the plot based on the
#' classifications. A suggested filtering threshold is provided for filtering
#' results based on the quality of the pattern-match.
#'
#' @param summaryList Classification summary table filtered to only include
#'   contigs with Prophage-like, Sloping and HighCovNoPattern classifications
#' @param suggFiltThresh TRUE or FALSE, Suggest a filtering threshold on the
#'   output pattern-match score histogram. Default is FALSE.
#' @return ggplot object
#' @keywords internal
resultsHisto <- function(summaryList, suggFiltThresh) {
    normMatchScore <- classifications <- NULL
    summaryTable <- summaryList[[1]]
    cleanSummaryTable <- summaryList[[2]]
    brks <- seq(
        min(cleanSummaryTable$normMatchScore),
        max(cleanSummaryTable$normMatchScore),
        length.out = 40
    )
    histogram <-
        hist(cleanSummaryTable$normMatchScore,
            breaks = brks,
            plot = FALSE
        )
    if (suggFiltThresh == TRUE) {
        max <- which.max(histogram$counts)
        SD <- sd(histogram$counts) * 2
        ST <- histogram$breaks[0:1 + (max + SD)][[1]]
        plot <- ggplot(data = cleanSummaryTable) +
            theme_bw() +
            geom_histogram(aes(x = normMatchScore, fill = classifications),
                breaks = brks
            ) +
            geom_vline(xintercept = ST, color = "red") +
            annotate(
                geom = "label",
                x = ST,
                y = max(histogram$counts),
                label = round(ST, digits = 2)
            ) +
            labs(
                title = "Distribution of Pattern-Match Classifications",
                x = "Normalized Pattern-Match Score",
                y = "count",
                caption = paste(
                    "(Lower scores are better pattern-matches) \n Suggested
                    Filtering Threshold=",
                    round(ST, digits = 2),
                    "\n",
                    (length(which(summaryTable[, 2] == "NoPattern"))),
                    "contigs were classified as 'NoPattern' and are not
                    displayed on this plot."
                )
            ) +
            theme(
                plot.title = element_text(size = 16),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                plot.caption = element_text(size = 12),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12)
            )
    } else {
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
    }
    return(plot)
}
