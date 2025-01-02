#' TrIdentSampleOutput
#'
#' The TrIdentClassifier output from the VLPFractionSamplePileup and
#' WholeCommunitySamplePileup files run with default parameters
#' Report...
#' @keywords internal
#' @usage data('TrIdentSampleOutput')
#' @format ## 'TrIdentSampleOutput' A list with 6 objects:
#' \describe{
#'    \item{SummaryTable}{A dataframe containing classifications for all
#'    contigs that were processed with pattern-matching}
#'    \item{CleanedSummaryTable}{SummaryTable dataframe filtered to remove
#'    contigs that recieved a 'None' classification}
#'    \item{PatternMatchInfo}{A list of lists containing pattern-match
#'    information for each classified contig}
#'    \item{FilteredOutContigTable}{A dataframe containing names of contigs
#'    that were filtered out prior to pattern-matching}
#'    \item{windowSize}{windowSize used in TrIdentClassifier function (1000)}
#'    \item{ResultHistogram}{a histogram displaying the overall abundance and
#'    quality of pattern-matches in addition to the composition of
#'    classifications. The displayed pattern-match scores are normalized by
#'    dividing each score by its associated contig length. The scores are
#'    normalized to visualize the overall quality of pattern-matching for the
#'    entire dataset.}
#'    }
#' @details A list object produced by the TrIdentClassifier function run on the
#' VLPFractionSamplePileup and WholeCommunitySamplePileup files run with default
#' parameters
#'
"TrIdentSampleOutput"
