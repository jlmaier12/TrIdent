% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patternCharacterizationFunctions.R
\name{prophageLikeElevation}
\alias{prophageLikeElevation}
\title{Determine Prophage-like read coverage elevation in whole-community}
\usage{
prophageLikeElevation(
  classifSummTable,
  prophageLikeClassifList,
  VLPpileup,
  WCpileup,
  windowSize,
  verbose
)
}
\arguments{
\item{classifSummTable}{Classification summary table}

\item{prophageLikeClassifList}{A list containing pattern match information
associated with all contigs classified as Prophage-like.}

\item{VLPpileup}{A table containing contig names, coverages averaged over 100
bp windows, and contig positions associated with mapping VLP-fraction reads
to whole-community contigs}

\item{WCpileup}{A table containing contig names, coverages averaged over 100
bp windows, and contig positions associated with mapping whole-community
reads to whole-community contigs}

\item{windowSize}{The window size used to re-average read coverage pileups}

\item{verbose}{TRUE or FALSE. Print progress messages to console. Default is
TRUE.}
}
\value{
dataframe
}
\description{
Determines whether a detected Prophage-like genetic element has read coverage
in the whole-community that is either elevated or depressed compared to the
average read coverage of the non-prophage region.
}
\keyword{internal}
