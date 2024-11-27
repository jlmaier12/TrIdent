#' Cleaned VLP-Fraction of Sample VLP-fraction pileup
#'
#' The VLPFractionSamplePileup after being 'cleaned' and reformatted for pattern-matching using the clean=TRUE parameter
#' Report...
#' @keywords internal
#' @usage data('CleanVLPFractionSamplePileup')
#' @format ## 'CleanVLPFractionSamplePileup'
#' A data frame with 25,668 rows and 3 columns:
#' \describe{
#'  \item{V1}{Contig accession}
#'  \item{V2}{Mapped read coverage averaged over a 100 bp window size}
#'  \item{V3}{Starting position (bp) of each 100 bp window. Restarts from 0 at the start of each new contig.}
#' }
#' @details
#' An example of a cleaned pileup file produced when clean=TRUE is used for
#' TrIdentClassifier, PlotTrIdentResults, and SpecializedTransductionID. While clean=TRUE is the default, if the user
#' selects clean=FALSE, their input pileup files must match the format of this example dataset.
#'
"CleanVLPFractionSamplePileup"
