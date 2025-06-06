% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fullSlopingPatternFunctions.R
\name{changeSlope}
\alias{changeSlope}
\title{Change slope of sloping pattern}
\usage{
changeSlope(
  leftOrRight,
  slopeBottom,
  halfToMaxReadCov,
  cov,
  viralSubset,
  windowSize
)
}
\arguments{
\item{leftOrRight}{Generate pattern for negative slope (left to right, i.e.
'Left') or positive slope (right to left, i.e. 'Right')}

\item{slopeBottom}{The value for the bottom of the sloping value}

\item{halfToMaxReadCov}{Half of the max VLP-fraction read coverage divided by
10}

\item{cov}{The value for the top of the slope}

\item{viralSubset}{A subset of the read coverage pileup that pertains only to
the contig currently being assessed}

\item{windowSize}{The window size used to re-average read coverage pileup}
}
\value{
List
}
\description{
Change the value of the slope used for the sloping pattern-match
}
\keyword{internal}
