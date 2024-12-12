#' Identify potential specialized transduction events on contigs classified as
#' Prophage-like
#'
#' Search contigs classified as Prophage-like for dense read coverage outside of
#' the pattern-match borders that may indicate specialized transduction. Returns
#' a list with the first object containing a summary table and the second object
#' containing a list of plots of with associated specialzied transduction search
#' results. If the plot is green, it has been identified as having potential
#' specialized transduction.
#'
#' @param VLPpileup VLP-fraction pileup file.
#' @param TrIdentResults Output from `TrIdentClassifier()`
#' @param noReadCov Number of basepairs of zero read coverage encountered before
#'   specialized transduction searching stops. Default is 500.
#' @param specTransLength Number of basepairs of non-zero read coverage needed
#'   for specialized transduction to be considered. Default is 2000.
#' @param specificContig Optional, Search a specific contig classified as
#'   Prophage-like ("NODE_1").
#' @param matchScoreFilter Optional, Filter plots using the normalized pattern
#'   match-scores. A suggested filtering threshold is provided by
#'   `TrIdentClassifier()` if `suggFiltThresh=TRUE`.
#' @param logScale TRUE or FALSE, display VLP-fraction read coverage in log10
#'   scale. Default is FALSE.
#' @param SaveFilesTo Provide a path to the directory you wish to save output
#'   to. `specializedTransductionID()` will make a folder within the provided
#'   directory to store results.
#' @return Large list containing two objects
#' @export
#' @examples
#' data("VLPFractionSamplePileup")
#' data("TrIdentSampleOutput")
#'
#' specTransduction <- specializedTransductionID(
#'     VLPpileup = VLPFractionSamplePileup,
#'     TrIdentResults = TrIdentSampleOutput
#' )
#'
#' specTransductionNODE62 <- specializedTransductionID(
#'     VLPpileup = VLPFractionSamplePileup,
#'     TrIdentResults = TrIdentSampleOutput,
#'     specificContig = "NODE_62"
#' )
specializedTransductionID <- function(VLPpileup,
                                      TrIdentResults,
                                      specificContig,
                                      noReadCov = 500,
                                      specTransLength = 2000,
                                      matchScoreFilter,
                                      logScale = FALSE,
                                      SaveFilesTo) {
    specTransInfo <- NULL
    TrIdentResultPatterns <- TrIdentResults[[3]]
    TrIdentResultSumm <- TrIdentResults[[1]]
    windowSize <- TrIdentResults[[5]]
    specTransSumm <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(specTransSumm) <-
        c(
            "contigName",
            "specTransduc",
            "left",
            "right",
            "lengthLeft",
            "lengthRight"
        )
    matchScoreFilter <-
        ifelse(missing(matchScoreFilter), Inf, matchScoreFilter)
    VLPpileup <- pileupFormatter(VLPpileup)
    specificContig <-
        ifelse(missing(specificContig), NA, specificContig)
    specTransCount <- 0
    plots <- list()
    J <- 1
    lapply(seq_along(TrIdentResultPatterns), function(i) {
        classification <- TrIdentResultPatterns[[i]][[7]]
        contigName <- TrIdentResultPatterns[[i]][[8]]
        normMatchScore <-
            TrIdentResultSumm[which(TrIdentResultSumm[, 1] ==
                contigName), 3]
        if (is.na(specificContig)) {
            if (classification == "Prophage-like" &
                normMatchScore < matchScoreFilter) {
                specTransInfo <<-
                  specTransductionSearch(
                        contigName,
                        VLPpileup,
                        TrIdentResultPatterns,
                        TrIdentResultSumm,
                        windowSize,
                        i,
                        noReadCov,
                        specTransLength,
                        logScale
                    )
                if (specTransInfo[[1]][[2]] == "yes") {
                    specTransCount <<- specTransCount + 1
                }
                specTransSumm[J, seq_len(6)] <<- specTransInfo[[1]]
                plots[[J]] <<- specTransInfo[[2]]
                J <<- J + 1
            }
        } else if (contigName == specificContig & classification ==
            "Prophage-like" &
            normMatchScore < matchScoreFilter) {
            specTransInfo <<- specTransductionSearch(
                contigName,
                VLPpileup,
                TrIdentResultPatterns,
                TrIdentResultSumm,
                windowSize,
                i,
                noReadCov,
                specTransLength,
                logScale
            )
            if (specTransInfo[[1]][[2]] == "yes") {
                specTransCount <<- specTransCount + 1
            }
            specTransSumm[J, seq_len(6)] <<- specTransInfo[[1]]
            plots[[J]] <<- specTransInfo[[2]]
            J <<- J + 1
        }
    })
    if (specTransCount == 0) {
        stop(
            "Selected contig is either not prophage-like, spelled incorrectly
            or has a match score below the chosen matchScoreFilter"
        )
    }
    message(
        specTransCount,
        " contigs have potential specialized transduction"
    )
    names(plots) <- specTransSumm[, 1]
    specTransList <- list(summaryTable = specTransSumm, Plots = plots)
    if (missing(SaveFilesTo) == FALSE) {
        ifelse(!dir.exists(paths = paste0(
            SaveFilesTo,
            "\\TrIdentSpecTransduction"
        )),
        dir.create(paste0(
            SaveFilesTo, "\\TrIdentSpecTransduction"
        )),
        stop(
            "'TrIdentSpecTransduction' folder exists already in the provided
            directory"
        )
        )
        lapply(
            names(plots),
            function(X) {
                ggsave(
                    filename = paste0(
                        SaveFilesTo,
                        "\\TrIdentSpecTransduction\\", X, ".png"
                    ),
                    plot = plots[[X]],
                    width = 8,
                    height = 4
                )
            }
        )
        write.table(
            specTransSumm,
            file = paste0(
                SaveFilesTo,
                "\\TrIdentSpecTransduction\\SpecTransducSummaryTable.csv"
            ),
            sep = ",",
            row.names = FALSE
        )
        return(specTransList)
    } else {
        return(specTransList)
    }
}
