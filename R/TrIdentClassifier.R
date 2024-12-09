#' Classify contigs as Prophage-like, Sloping,
#' HighCovNoPattern, and NoPattern
#'
#' Performs all the pattern-matching and summarizes the
#' results into a list. The first item in the list is a table
#' consisting of the summary information of all the contigs that
#' passed through pattern-matching (i.e were not filtered out).
#' The second item in the list is a table consisting of the
#' summary information of all contigs that were classified via
#' pattern-matching. The third item in the list contains the
#' pattern-match information associated with each contig in the
#' previous table. The fourth object in the list is a table containing
#' the contigs that were filtered out prior to pattern-matching. The
#' fifth item is the windowSize used for the search.
#'
#' @param VLPpileup VLP-fraction pileup file.
#' @param WCpileup A whole-community pileup file.
#' @param windowSize
#'  The number of basepairs to average read coverage
#'  values over. Options are 100,
#'  200, 500, 1000 ONLY.  Default is 1000.
#' @param minBlockSize
#'  The minimum size (in bp) of the Prophage-like block
#'  pattern. Default is 10000.
#' @param maxBlockSize
#'  The maximum size (in bp) of the Prophage-like block
#'  pattern. Default is NA
#'  (no maximum).
#' @param minContigLength
#'  The minimum contig size (in bp) to perform
#'  pattern-matching on. Must be at
#'  least 25000. Default is 30000.
#' @param minSlope
#'  The minimum slope value to test for sloping patterns.
#'  Default is 0.001 (i.e minimum change of 10x read
#'  coverage over 100,000 bp).
#' @param SaveFilesTo
#'  Optional, Provide a path to the directory you wish to
#'  save output to. A folder will be made within the provided
#'  directory to store results.
#' @param suggFiltThresh
#'  TRUE or FALSE, Suggest a filtering threshold for TrIdent
#'  classifications based on the normalized pattern-match
#'  scores. Default is FALSE.
#' @importFrom utils capture.output
#' @return Large list containing 5 objects
#' @export
#'
#' @examples
#' data("VLPFractionSamplePileup")
#' data("WholeCommunitySamplePileup")
#' TrIdent_results <- TrIdentClassifier(
#'     VLPpileup = VLPFractionSamplePileup,
#'     WCpileup = WholeCommunitySamplePileup
#' )
TrIdentClassifier <- function(VLPpileup,
                                WCpileup,
                                windowSize = 1000,
                                minBlockSize = 10000,
                                maxBlockSize = Inf,
                                minContigLength = 30000,
                                minSlope = 0.001,
                                suggFiltThresh = FALSE,
                                SaveFilesTo) {
    ## error catching
    if (!(windowSize %in% list(100, 200, 500, 1000))) {
        stop("windowSize must be either 100, 200, 500, or 1000 bp!")
    }
    if (minContigLength < 25000) {
        stop("minContigLength must be at least 25,000 bp for pattern-matching!")
    }
    if (minBlockSize <= 1000) {
        stop("minBlockSize must be greater than 1000 bp!")
    }
    if (nrow(VLPpileup) != nrow(WCpileup)) {
        stop("VLP and WC pileup files have differing row numbers")
    }
    if (abs(VLPpileup[1, 3] - VLPpileup[2, 3]) != 100 |
        abs(WCpileup[1, 3] - WCpileup[2, 3]) != 100) {
        stop("pileup files MUST have a windowSize/binsize of 100!")
    }

    ## main algorithm start
    startTime <- Sys.time()
    message("Reformatting pileup files")
    VLPpileup <- pileupFormatter(VLPpileup)
    WCpileup <- pileupFormatter(WCpileup)
    message("Starting pattern-matching...")
    classificationSummary <-
        patternMatcher(
            VLPpileup,
            WCpileup,
            windowSize,
            minBlockSize,
            maxBlockSize,
            minContigLength,
            minSlope
        )
    summaryTable <- contigClassSumm(classificationSummary[[1]])
    summaryTable <- VLPtoWCRatioCalc(summaryTable, WCpileup, VLPpileup)
    summaryTable <-
        patternMatchSize(
            summaryTable, classificationSummary[[1]],
            windowSize
        )
    summaryTable <- prophageLikeElevation(
        summaryTable,
        allProphageLikeClassifs(classificationSummary[[1]]),
        VLPpileup,
        WCpileup,
        windowSize
    )
    summaryTable <- slopeSumm(
        summaryTable,
        allSlopingClassifs(classificationSummary[[1]])
    )

    message("Finalizing output")
    summaryList <- list(
        SummaryTable = summaryTable,
        CleanedSummaryTable = summaryTable[which(
            summaryTable[, 2] ==
                "Prophage-like" |
                summaryTable[, 2] ==
                    "Sloping" |
                summaryTable[, 2] ==
                    "HighCovNoPattern"
        ), ],
        PatternMatchInfo = allPatternMatches(
            classificationSummary[[1]],
            summaryTable
        ),
        FilteredOutContigTable = classificationSummary[[2]],
        windowSize = windowSize
    )
    endTime <- Sys.time()
    duration <- difftime(endTime, startTime)
    message("Execution time: ", round(duration[[1]], 2), units(duration))
    message(
        length(which(
            classificationSummary[[2]][, 2] == "Low VLP-fraction read cov"
        )),
        " contigs were filtered out based on low read coverage"
    )
    message(
        length(which(
            classificationSummary[[2]][, 2] == "Contig length too small"
        )),
        " contigs were filtered out based on length"
    )
    table <- (table(summaryList[[1]][, 2]))
    message(paste0(capture.output(table), collapse = "\n"))
    message(
        length(which(summaryList[[1]][, 8] == "Elevated")),
        " of the prophage-like classifications are highly active or abundant"
    )
    message(
        length(which(summaryList[[1]][, 8] == "Depressed")),
        " of the prophage-like classifications are mixed, i.e. heterogenously
        integrated into their bacterial host population"
    )
    plot(resultsHisto(summaryList, suggFiltThresh))
    if (missing(SaveFilesTo) == FALSE) {
        ifelse(!dir.exists(paths = paste0(SaveFilesTo, "\\TrIdentOutput")),
            dir.create(paste0(SaveFilesTo, "\\TrIdentOutput")),
            stop(
                "'TrIdentOutput' folder exists already in the provided directory"
            )
        )
        write.table(
            summaryTable,
            file = paste0(
                SaveFilesTo,
                "\\TrIdentOutput\\TrIdentSummaryTable.csv"
            ),
            sep = ",",
            row.names = FALSE
        )
        write.table(
            summaryTable[which(
                summaryTable[, 2] ==
                    "Prophage-like" |
                    summaryTable[, 2] ==
                        "Sloping" | summaryTable[, 2] ==
                    "HighCovNoPattern"
            ), ],
            file = paste0(
                SaveFilesTo,
                "\\TrIdentOutput\\TrIdentSummaryTableCleaned.csv"
            ),
            sep = ",",
            row.names = FALSE
        )
        write.table(
            classificationSummary[[2]],
            file = paste0(
                SaveFilesTo,
                "\\TrIdentOutput\\TrIdentFilteredOutContigs.csv"
            ),
            sep = ",",
            row.names = FALSE
        )
        ggsave(
            filename = paste0(
                SaveFilesTo,
                "\\TrIdentOutput\\MatchScoreDensityPlot.png"
            ),
            plot = plot,
            width = 4,
            height = 4
        )
        return(summaryList)
    } else {
        return(summaryList)
    }
}
