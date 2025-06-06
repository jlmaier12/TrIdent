% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UtilFunctions.R
\name{pileupFormatter}
\alias{pileupFormatter}
\title{Correctly formats pileup files.}
\usage{
pileupFormatter(pileup)
}
\arguments{
\item{pileup}{A table containing contig names, read coverages averaged over
100 bp windows,and contig positions}
}
\value{
dataframe
}
\description{
Places columns in correct order and renames columns. Cleans the contig labels
to remove excess information after whitespace.
}
\keyword{internal}
