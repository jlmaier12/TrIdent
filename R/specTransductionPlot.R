#' Specialized transduction plot
#'
#' Plot search results of `specializedTransductionID()`
#'
#' @param viralSubsetZoom contig subset surrounding Prophage-like pattern-match
#' @param startPosBp Left border position
#' @param endPosBp Right border position
#' @param SpecTransLeft End position of spec transduction on left border
#' @param specTransRight End position of spec transduction on right border
#' @param contigName The reference name of the contig currently being assessed
#'   (i.e "NODE_1")
#' @param classifPatternMatches The pattern match information associated with
#'   each contig classified as prophage-like,  sloping, or HighCovNoPattern
#' @param i The index for the contig currently being assessed
#' @param specTransSumm Results for spec transduction search
#' @param logScale If TRUE, coverage is plotted in log10. If FALSE, raw coverage
#'   values are plotted. Default is FALSE.
#' @param classifSumm The summary information associated with each contig
#'   classified as Prophage-like, Sloping, or HighCovNoPattern
#' @return ggplot object
#' @keywords internal
specTransductionPlot <-
    function(viralSubsetZoom,
             startPosBp,
             endPosBp,
             SpecTransLeft,
             specTransRight,
             contigName,
             classifPatternMatches,
             i,
             specTransSumm,
             logScale,
             classifSumm) {
        position <- logcoverage <- NULL
        if (classifSumm[
            which(classifSumm[, 1] == contigName),
            6
        ] == "Elevated") {
            prophageLikeInfo <- "Highly active/abundant prophage-like element"
        } else if (classifSumm[
            which(classifSumm[, 1] == contigName),
            6
        ] == "Depressed") {
            prophageLikeInfo <-
                "Not homogenously present/integrated prophage-like element"
        } else {
            prophageLikeInfo <- NULL
        }
        fill <-
            ifelse(specTransSumm[2] == "yes", "seagreen", "deepskyblue3")
        alphaL <- ifelse(specTransSumm[3] == "yes", 1, 0)
        alphaR <- ifelse(specTransSumm[4] == "yes", 1, 0)
        coverageType <-
            if (logScale == TRUE) {
                viralSubsetZoom$logcoverage
            } else {
                viralSubsetZoom$coverage
            }
        plot <- (
            ggplot(data = as.data.frame(viralSubsetZoom), aes(
                x = position,
                y = coverageType
            )) +
                geom_area(fill = fill) +
                geom_vline(
                    xintercept = c(startPosBp, endPosBp),
                    linewidth = 1
                ) +
                geom_vline(
                    xintercept = SpecTransLeft,
                    color = "red",
                    alpha = alphaL,
                    linewidth = 1
                ) +
                geom_vline(
                    xintercept = specTransRight,
                    color = "red",
                    alpha = alphaR,
                    linewidth = 1
                ) +
                scale_x_continuous(expand = c(0, 0)) +
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.subtitle = element_text(size = 10),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    text = element_text(size = 15)
                ) +
                labs(
                    title = paste(
                        contigName, classifPatternMatches[[i]][[8]],
                        prophageLikeInfo
                    ),
                    subtitle = paste(
                        "Specialized transduction on left:",
                        specTransSumm[3],
                        ",",
                        "on right:",
                        specTransSumm[4]
                    ),
                    x = "Contig Position (bp)",
                    y = paste(
                        "VLP-Fraction Read Coverage",
                        ifelse(logScale == TRUE,
                            "(Log10)", ""
                        )
                    )
                )
        )
        return(plot)
    }
