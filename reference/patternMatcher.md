# Main pattern-matching function

Creates the viralSubset, representative of one contig, that is used as
input for each individual pattern-matching function. After the
information associated with the best match for each pattern is obtained,
the pattern with the smallest match score is used to classify the contig
being assessed. Prior to the pattern-matching, contigs smaller than the
minContigLength and contigs without 5,000 bp of 10x read coverage are
removed.

## Usage

``` r
patternMatcher(
  VLPpileup,
  WCpileup,
  windowSize,
  minBlockSize,
  maxBlockSize,
  minContigLength,
  minSlope,
  verbose
)
```

## Arguments

- VLPpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping VLP-fraction
  reads to whole-community contigs

- WCpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping whole-community
  reads to whole-community contigs

- windowSize:

  The window size used to re-average read coverage datasets

- minBlockSize:

  The minimum size of the prophage-like block pattern. Default is 10,000
  bp.

- maxBlockSize:

  The maximum size of the prophage-like block pattern. Default is NA

- minContigLength:

  The minimum contig size (in bp) to perform pattern-matching on. Must
  be at least 20,000 bp. Default is 30,000 bp.

- minSlope:

  The minimum slope value to test for sloping patterns

- verbose:

  TRUE or FALSE. Print progress messages to console. Default is TRUE.

## Value

List containing three objects.
