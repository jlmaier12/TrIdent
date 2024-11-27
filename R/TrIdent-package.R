#' @title \bold{TrIdent} - \bold{Tr}ansduction \bold{Ident}ification
#'
#' @description
#' Automatic detection, classification and characterization of transduction
#' events in transductomics datasets using read coverage pattern-matching.
#'
#' Please see \href{https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5}{Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities}
#' for more information on the transductomics method, data and analysis workflow.
#'
#' @details
#' The three main functions in TrIdent are:
#' \enumerate{
#' \item \code{\link{TrIdentClassifier}} performs the pattern-matching, classification and
#' characterization of read coverage patterns on contigs.
#' \item \code{\link{plotTrIdentResults}} plots the results from \code{TrIdentClassifier()}
#' \item \code{\link{specializedTransductionID}} searches contigs classified as Prophage-like
#' by \code{TrIdentClassifier()} for potential specialized transduction
#'}
#'
#' @author Jessie Maier \email{jlmaier@ncsu.edu} & Jorden Rabasco \email{jrabasc@ncsu.edu}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @import stringr
#' @import tidyr
#' @importFrom graphics hist
#' @importFrom roll roll_sd
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom utils capture.output
#' @importFrom utils write.table
## usethis namespace: end
NULL
