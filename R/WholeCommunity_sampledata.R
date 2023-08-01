#' Whole-Community Fraction of Sample Dataset
#'
#' A subset of contigs from the raw whole-community fraction read coverage pileup file generated from BBMap's pileup.sh.
#' Report...
#'
#' @format ## 'WholeCommunity_sampledata'
#' A data frame with 106,167 rows and 4 columns:
#' \describe{
#'  \item{V1}{Contig accession}
#'  \item{V2}{Read coverage binned over binsize used by pileup.sh}
#'  \item{V3}{Starting position (bp) of each bin used by pileup.sh. Restarts from 0 at the start of each new contig.}
#'  \item{V4}{Starting position (bp) of each bin used by pileup.sh. Does NOT restart at the start of each new contig.}
#' }

"WholeCommunity_sampledata"
