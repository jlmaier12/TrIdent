#' Classify contigs as Prophage-like, Sloping, HighCovNoPattern, and NoPattern
#'
#' Performs all the pattern-matching and summarizes the results into a list. The
#' first item in the list is a table consisting of the summary information of
#' all the contigs that passed through pattern-matching (i.e were not filtered
#' out). The second item in the list is a table consisting of the summary
#' information of all contigs that were classified via pattern-matching. The
#' third item in the list contains the pattern-match information associated with
#' each contig in the previous table. The fourth object in the list is a table
#' containing the contigs that were filtered out prior to pattern-matching. The
#' fifth item is the windowSize used for the search.
#'
#' @param VLPpileup VLP-fraction pileup file generated by mapping sequencing reads 
#' from a sample's ultra-purified VLP-fraction mapped to the sample's whole-community
#' metagenome assembly. The pileup file MUST have the following format:
#' * V1: Contig accession
#' * V2: Mapped read coverage values averaged over 100 bp windows
#' * V3: Starting position (bp) of each 100 bp window. Restarts from 0 at the
#' start of each new contig.
#' * V4: Starting position (bp) of each 100 bp window. Does NOT restart at the
#' start of each new contig.
#' @param WCpileup A whole-community pileup file generated by mapping sequencing reads
#' from a sample's whole-community mapped to the sample's whole-community metagenome assembly.
#' The pileup file MUST have the following format:
#' * V1: Contig accession
#' * V2: Mapped read coverage values averaged over 100 bp windows
#' * V3: Starting position (bp) of each 100 bp window. Restarts from 0 at the
#' start of each new contig.
#' * V4: Starting position (bp) of each 100 bp window. Does NOT restart at the
#' start of each new contig.
#' @param windowSize The number of basepairs to average read coverage values
#'   over. Options are 100, 200, 500, 1000 ONLY.  Default is 1000.
#' @param minBlockSize The minimum size (in bp) of the Prophage-like block
#'   pattern. Default is 10000. Must be at least 1000.
#' @param maxBlockSize The maximum size (in bp) of the Prophage-like block
#'   pattern. Default is NA (no maximum).
#' @param minContigLength The minimum contig size (in bp) to perform
#'   pattern-matching on. Must be at least 25000. Default is 30000.
#' @param minSlope The minimum slope value to test for sloping patterns. Default
#'   is 0.001 (i.e minimum change of 10x read coverage over 100,000 bp).
#' @param VLPReads Optional, the number of VLP-fraction reads used for mapping 
#'   and creation of pileup.
#' @param WCReads Optional, the number of WC reads used for mapping and 
#'   creation of pileup.
#' @param SaveFilesTo Optional, Provide a path to the directory you wish to save
#'   output to. A folder will be made within the provided directory to store
#'   results.
#' @param verbose TRUE or FALSE. Print progress messages to console. Default is TRUE.
#' @importFrom utils capture.output
#' @return Large list containing 5 objects
#' @export
#' @examples
#' data("VLPFractionSamplePileup")
#' data("WholeCommunitySamplePileup")
#'
#' TrIdent_results <- TrIdentClassifier(
#'   VLPpileup = VLPFractionSamplePileup,
#'   WCpileup = WholeCommunitySamplePileup
#' )
TrIdentClassifier <- function(VLPpileup,
                              WCpileup,
                              windowSize = 1000,
                              minBlockSize = 10000,
                              maxBlockSize = Inf,
                              minContigLength = 30000,
                              minSlope = 0.001,
                              VLPReads,
                              WCReads,
                              verbose = TRUE,
                              SaveFilesTo) {
  ## error catching
  if (!(windowSize %in% list(100, 200, 500, 1000))) {
    stop("windowSize must be either 100, 200, 500, or 1000 bp!")
  }
  if (minContigLength <= 25000) {
    stop("minContigLength must be at least 25,000 bp for pattern-matching!")
  }
  if (minBlockSize <= 1000) {
    stop("minBlockSize must be greater than 1000 bp!")
  }
  ## input validation
  if (nrow(VLPpileup) != nrow(WCpileup)) {
    stop("VLP and WC pileup files have differing row numbers")
  }
  if (abs(VLPpileup[1, 3] - VLPpileup[2, 3]) != 100 |
    abs(WCpileup[1, 3] - WCpileup[2, 3]) != 100) {
    stop("pileup files MUST have a windowSize/binsize of 100!")
  }
  if (all(VLPpileup[, 1] == WCpileup[, 1]) == FALSE) {
    stop("The first column of the VLP and WC pileup file should be identical if
         mapping was performed correctly...")
  }
  if (missing(VLPReads) & missing(WCReads) == FALSE |
      missing(VLPReads) == FALSE & missing(WCReads)) {
      stop("Both the WC and VLP-fraction read coverage values must be supplied.
           You can not supply one without the other.")
  }
  ## main algorithm start
  startTime <- Sys.time()
  if (verbose) {message("Reformatting pileup files")}
  VLPpileup <- pileupFormatter(VLPpileup)
  WCpileup <- pileupFormatter(WCpileup)
  if (verbose) {message("Starting pattern-matching...")}
  classificationSummary <-
    patternMatcher(
      VLPpileup,
      WCpileup,
      windowSize,
      minBlockSize,
      maxBlockSize,
      minContigLength,
      minSlope,
      verbose
    )
  VLPReads <- ifelse(missing(VLPReads), 1, VLPReads)
  WCReads <- ifelse(missing(WCReads), 1, WCReads)
  summaryTable <- contigClassSumm(classificationSummary[[1]])
  summaryTable <- VLPtoWCRatioCalc(summaryTable, WCpileup, VLPpileup, VLPReads, WCReads)
  summaryTable <-
    patternMatchSize(
      summaryTable, classificationSummary[[1]],
      windowSize,
      verbose
    )
  summaryTable <- prophageLikeElevation(
    summaryTable,
    allProphageLikeClassifs(classificationSummary[[1]]),
    VLPpileup,
    WCpileup,
    windowSize,
    verbose
  )
  summaryTable <- slopeSumm(
    summaryTable,
    allSlopingClassifs(classificationSummary[[1]]),
    windowSize
  )
  if (verbose) {message("Finalizing output")}
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
  plot <- ifelse(nrow(summaryList[[2]]) > 0, resultsHisto(summaryList), NULL)
  summaryList <- c(summaryList, ResultHistogram = plot)
  endTime <- Sys.time()
  duration <- difftime(endTime, startTime)
  if (verbose) {
      message(
          "Execution time: ", round(duration[[1]], 2), units(duration))}
  if (verbose) {
      message(
          length(which(
              classificationSummary[[2]][, 2] == "Low VLP-fraction read cov"
      )),
      " contigs were filtered out based on low read coverage"
    )
  }
  if (verbose) {
    message(
      length(which(
        classificationSummary[[2]][, 2] == "Contig length too small"
      )),
      " contigs were filtered out based on length"
    )
  }
  table <- (table(summaryList[[1]][, 2]))
  if (verbose) {
    message(paste0(capture.output(table), collapse = "\n"))
  }
  if (verbose) {
    message(
      length(which(summaryList[[1]][, 8] == "Elevated")),
      " of the prophage-like classifications are highly active or abundant"
    )
  }
  if (verbose) {
    message(
      length(which(summaryList[[1]][, 8] == "Depressed")),
      " of the prophage-like classifications are mixed, i.e. heterogenously
        integrated into their bacterial host population"
    )
  }
  if (missing(SaveFilesTo) == FALSE) {
    ifelse(!dir.exists(paths = paste0(SaveFilesTo, "\\TrIdentOutput")),
      dir.create(paste0(SaveFilesTo, "\\TrIdentOutput")),
      stop(
        "'TrIdentOutput' folder exists already
                in the provided directory"
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
        "\\TrIdentOutput\\TrIdentResultHistogram.png"
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
