# Sloping pattern builder

Build a sloping pattern that consists of a sloping line spanning the
contig being assessed. The line slopes from left to right. The slope of
the line is changed, but the pattern is not translated across the
contig.

## Usage

``` r
fullSlope(
  viralSubset,
  windowSize,
  minSlope,
  minSlopeSize,
  searchMethod,
  DirectMaxEval,
  globalLocal
)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileup

- minSlope:

  The minimum slope value to test for sloping patterns

- minSlopeSize:

  The minimum width of sloping patterns.

- searchMethod:

  Search method to use. Either "grid" for the original grid search or
  "direct" for DIRECT global optimization.

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make.

- globalLocal:

  Use global or local DIRECT search. Default is local.

## Value

List containing two objects
