% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fullSlopingPatternFunctions.R
\name{fullSlope}
\alias{fullSlope}
\title{Sloping pattern builder}
\usage{
fullSlope(viralSubset, windowSize, minSlope)
}
\arguments{
\item{viralSubset}{A subset of the read coverage pileup that pertains only to
the contig currently being assessed}

\item{windowSize}{The window size used to re-average read coverage pileup}

\item{minSlope}{The minimum slope value to test for sloping patterns}
}
\value{
List containing two objects
}
\description{
Build a sloping pattern that consists of a sloping line spanning the contig
being assessed. The line slopes from left to right. The slope of the line is
changed, but the pattern is not translated across the contig.
}
\keyword{internal}
