% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VLPFractionSamplePileup.R
\docType{data}
\name{VLPFractionSamplePileup}
\alias{VLPFractionSamplePileup}
\title{VLP-Fraction of Sample Dataset}
\format{
## 'VLPFractionSamplePileup' A data frame with 10,805 rows and 4
  columns:
\describe{
 \item{V1}{Contig accession}
 \item{V2}{Mapped read coverage averaged over a 100 bp window size}
 \item{V3}{Starting position (bp) of each 100 bp window. Restarts from 0 at
 the start of each new contig.}
 \item{V4}{Starting position (bp) of each 100 bp window. Does NOT restart at
 the start of each new contig.}
}
}
\source{
<https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5>
}
\usage{
data('VLPFractionSamplePileup')
}
\description{
A subset of contigs from the raw VLP-fraction read coverage pileup file
generated from BBMap's pileup.sh. Report...
}
\details{
This dataset represents one half of a complete transductomics
dataset which is comprised of two parts-a whole-community fraction and a
viral-like particle (VLP)-fraction. This dataset represents the VLP fraction
and was generated by purifying VLPs from a conventional mouse fecal
homogenate using CsCl density gradient ultracentrifugation. The VLP-fraction
extracted DNA was sequenced with Illumina (paired-end mode, 150 bp reads) The
sequencing reads were mapped to the associated whole-community assembly using
BBMap. The bbmap.sh bincov parameter with covbinsize=100 was used to create a
pileup file with 100 bp windows. A subset of 10 contigs from the pileup file
were selected for this sample dataset. The contigs were chosen because their
associated read coverage patterns in the VLP-fraction exemplify TrIdent's
pattern-matching functionality across classifications:
NODE_617:Prophage-like, active/abundant, with spec transduction
NODE_135:Prophage-like, off one side of contig, no spec transduction
NODE_352:Sloping, left to right slope
NODE_1088: Sloping, right to left slope
NODE_2060: Sloping, right to left slope with start
NODE_1401: None, no pattern match
NODE_62: Prophage-like, with spec transduction
NODE_368: Prophage-like, not homogeneously integrated/present, no spec
transduction
NODE_560: HighCovNoPattern
NODE_1165: None, filtered out
To access the sequencing data used to generate this pileup file and for
additional details on the assembly and mapping parameters, refer to the
reference below:
Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al. Transductomics:
sequencing-based detection and analysis of transduced DNA in pure cultures
and microbial communities. Microbiome 8, 158 (2020).
https://doi.org/10.1186/s40168-020-00935-5
}
\keyword{internal}
