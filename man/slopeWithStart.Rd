% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partialSlopingPatternFunctions.R
\name{slopeWithStart}
\alias{slopeWithStart}
\title{Sloping pattern with an initial jump-up in read coverage}
\usage{
slopeWithStart(viralSubset, windowSize, minSlope)
}
\arguments{
\item{viralSubset}{A subset of the read coverage pileup that pertains only to
the contig currently being assessed}

\item{windowSize}{The window size used to re-average read coverage pileups}

\item{minSlope}{The minimum slope value to test for sloping patterns}
}
\value{
List containing two objects
}
\description{
Build, translate, and change slope of sloping pattern with slope start
}
\keyword{internal}
