#' Whole-Community Fraction of Sample Dataset
#'
#' A subset of contigs from the raw whole-community fraction read coverage pileup file generated during read mapping.
#' Report...
#'
#' @keywords internal
#' @usage data('WholeCommunitySamplePileup')
#' @format ## 'WholeCommunitySamplePileup'
#' A data frame with 10,805 rows and 4 columns:
#' \describe{
#'  \item{V1}{Contig accession}
#'  \item{V2}{Mapped read coverage averaged over a 100 bp window size}
#'  \item{V3}{Starting position (bp) of each 100 bp window. Restarts from 0 at the start of each new contig.}
#'  \item{V4}{Starting position (bp) of each 100 bp window. Does NOT restart at the start of each new contig.}
#' }
#' @details
#' This dataset represents one half of a complete transductomics dataset which is comprised of two parts-
#' a whole-community fraction and a viral-like particle (VLP)-fraction. This dataset represents the
#' whole-community fraction and was generated from a conventional mouse fecal homogenate. The whole-community
#' extracted DNA was sequenced with Illumina (paired-end mode, 150 bp reads)
#' after which the metagenome was assembled. The sequencing reads were mapped to the assembled contigs
#' using BBMap from the BBTools suite. The bbmap.sh bincov parameter with covbinsize=100 was used to create a pileup file with 100 bp windows. A
#' subset of 10 contigs from the pileup file were selected for this sample dataset. The contigs were chosen because their associated
#' read coverage patterns in the VLP-fraction exemplify TrIdent's pattern-matching and characterization functionality across classifications:
#' NODE_617: Prophage-like, active/abundant, with spec transduction
#' NODE_125: Prophage-like, off one side of contig, no spec transduction
#' NODE_352: Sloping, left to right slope
#' NODE_251: Sloping, right to left slope
#' NODE_2060: Sloping, right to left slope with start
#' NODE_1401: None, no pattern match
#' NODE_62: Prophage-like, with spec transduction
#' NODE_368: Prophage-like, not homogeneously integrated/present, no spec transduction
#' NODE_560: HighCovNoPattern
#' NODE_1165: None, filtered out
#' To access the sequencing data used to generate this pileup file and for additional details
#' on the assembly and mapping parameters, refer to the reference below:
#' Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al.
#' Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities.
#' Microbiome 8, 158 (2020). https://doi.org/10.1186/s40168-020-00935-5
#' @source <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5>
"WholeCommunitySamplePileup"
