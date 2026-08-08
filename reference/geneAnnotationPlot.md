# Gene annotation plot

Plot read coverage and location of gene annotations that match the
keywords and search criteria for contig currently being assessed

## Usage

``` r
geneAnnotationPlot(
  geneAnnotSubset,
  keywords,
  VLPpileupSubset,
  colIdx,
  startbpRange,
  endbpRange,
  pattern,
  windowSize
)
```

## Arguments

- geneAnnotSubset:

  Subset of gene annotations to be plotted

- keywords:

  The key-word(s) used for the search.

- VLPpileupSubset:

  A subset of the pileup associated with the contig/chunk being assessed

- colIdx:

  The gff column being searched

- startbpRange:

  The basepair at which the search is started if a 'specific' search is
  used

- endbpRange:

  The basepair at which the search is ended if a 'specific' search is
  used

- pattern:

  The pattern-match information associated with the contig/chunk being
  assessed

- windowSize:

  The number of basepairs to average read coverage values over.
