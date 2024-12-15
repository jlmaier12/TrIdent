#' Plot read coverage graphs of contigs classified as Prophage-like, Sloping, or
#' HighCovNoPattern
#'
#' Plot the read coverages of a contig and its associated pattern-match for
#' Prophage-like, Sloping and HighCovNoPattern classifications. Returns a list
#' of ggplot objects.
#'
#' @param VLPpileup VLP-fraction pileup file.
#' @param WCpileup Whole-community pileup file.
#' @param TrIdentResults Output from `TrIdentClassifier()`.
#' @param matchScoreFilter Optional, Filter plots using the normalized pattern
#'   match-scores. A suggested filtering threshold is provided by
#'   `TrIdentClassifier()` if `suggFiltThresh=TRUE`.
#' @param saveFilesTo Optional, Provide a path to the directory you wish to save
#'   output to. A folder will be made within the provided directory to store
#'   results.
#' @return Large list containing ggplot objects
#' @export
#' @examples
#' data("VLPFractionSamplePileup")
#' data("WholeCommunitySamplePileup")
#' data("TrIdentSampleOutput")
#' patternMatches <- plotTrIdentResults(
#'     VLPpileup = VLPFractionSamplePileup,
#'     WCpileup = WholeCommunitySamplePileup,
#'     TrIdentResults = TrIdentSampleOutput
#' )
plotTrIdentResults <- function(VLPpileup,
                               WCpileup,
                               TrIdentResults,
                               matchScoreFilter,
                               saveFilesTo) {
    position <- coverage <- NULL
    windowSize <- TrIdentResults[[5]]
    cleanSummaryTable <- TrIdentResults[[3]]
    summaryTable <- TrIdentResults[[1]]
    MSF <-
        ifelse(missing(matchScoreFilter) == TRUE, 0, matchScoreFilter)
    VLPpileup <- pileupFormatter(VLPpileup)
    WCpileup <- pileupFormatter(WCpileup)
    plots <- list()
    contigNames <- c()
    plots <- lapply(seq_along(cleanSummaryTable), function(i) {
        contigName <- cleanSummaryTable[[i]][[8]]
        viralSubset <-
            changeWindowSize(
                VLPpileup[which(VLPpileup[, 1] == contigName), ],
                windowSize
            )
        microbialSubset <-
            changeWindowSize(
                WCpileup[which(WCpileup[, 1] == contigName), ],
                windowSize
            )
        patternMatchInfo <-
            summaryTable[which(summaryTable[, 1] == contigName), ]
        classification <- patternMatchInfo[, 2]
        pattern <-
            patternBuilder(viralSubset, cleanSummaryTable, classification, i)
        patternMatch <- cbind(viralSubset, pattern)
        matchLength <- patternMatchInfo[, 5]
        matchscoreQC <-
            (cleanSummaryTable[[i]][[1]]) / mean(viralSubset$coverage)
        if (MSF != 0) {
            if (matchscoreQC > MSF) {
                return()
            }
        }
        if (classification == "Sloping") {
            subtitleInfo <- paste(
                "Slope:",
                patternMatchInfo[10]
            )
        } else if (classification == "HighCovNoPattern") {
            subtitleInfo <- paste("VLP:WC ratio:", patternMatchInfo[4])
        } else if (classification == "Prophage-like") {
            if (is.na(patternMatchInfo[8]) == TRUE) {
                subtitleInfo <- NULL
            } else if (patternMatchInfo[8] == "Elevated") {
                subtitleInfo <- "Active/highly abundant Prophage-like element"
            } else if (patternMatchInfo[8] == "Depressed") {
                subtitleInfo <- "Not homogenously integrated Prophage-like
                element"
            } else {
                subtitleInfo <- NULL
            }
        }
        wholecomm_plot <-
            ggplot(data = microbialSubset, aes(x = position, y = coverage)) +
            geom_area(fill = "deepskyblue3") +
            labs(
                title = paste(contigName, "Classification:", classification),
                subtitle = paste(
                    "Match size (bp):", matchLength,
                    subtitleInfo
                ),
                x = " ",
                y = "Whole-community \n read coverage"
            ) +
            scale_x_continuous(expand = c(0, 0)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 14),
                axis.text = element_text(size = 11),
                plot.subtitle = element_text(size = 12),
                plot.title = element_text(size = 14),
                plot.margin = margin(
                    t = 0,
                    r = 6,
                    b = 0,
                    l = 2
                )
            )
        Overlay_plot <-
            ggplot(data = patternMatch, aes(x = position, y = coverage)) +
            geom_area(fill = "deepskyblue3") +
            geom_line(aes(y = pattern),
                color = "black",
                linewidth = 1
            ) +
            labs(
                x = "Contig position (bp)",
                y = "VLP-fraction \n read coverage"
            ) +
            scale_x_continuous(expand = c(0, 0)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 14),
                axis.text = element_text(size = 11),
                plot.margin = margin(
                    t = 0,
                    r = 6,
                    b = 0,
                    l = 2
                )
            )
        contigNames <<- c(contigNames, contigName)
        combined_plot <- (wholecomm_plot / Overlay_plot)
        combined_plot
    })
    plots <- Filter(Negate(is.null), plots)
    names(plots) <- contigNames
    if (missing(saveFilesTo) == FALSE) {
        ifelse(!dir.exists(paths = paste0(saveFilesTo, "\\TrIdentOutput")),
            dir.create(paste0(saveFilesTo, "\\TrIdentOutput")),
            stop(
                "'TrIdentOutput' folder exists already in the provided
                directory"
            )
        )
        lapply(
            names(plots),
            function(X) {
                ggsave(
                    filename = paste0(
                        saveFilesTo,
                        "\\TrIdentOutput\\", X, ".png"
                    ),
                    plot = plots[[X]],
                    width = 8,
                    height = 4
                )
            }
        )
        return(plots)
    } else {
        return(plots)
    }
}
