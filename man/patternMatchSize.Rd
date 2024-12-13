% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patternMatchSize.R
\name{patternMatchSize}
\alias{patternMatchSize}
\title{Pattern-match size calculator}
\usage{
patternMatchSize(classifSumm, classifList, windowSize)
}
\arguments{
\item{classifSumm}{Classification summary table}

\item{classifList}{A list containing pattern match information associated
with all contig classifications}

\item{windowSize}{The window size used to re-average read coverage pileups}
}
\value{
dataframe
}
\description{
Calculate the size (bp) of the matching region for Prophage-like and Sloping
patterns
}
\keyword{internal}
