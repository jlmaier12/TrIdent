# Sloping pattern builder

Build a sloping pattern that consists of a sloping line spanning the
contig being assessed. The line slopes from left to right. The slope of
the line is changed, but the pattern is not translated across the
contig.

## Usage

``` r
fullSlope(viralSubset, windowSize, minSlope)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileup

- minSlope:

  The minimum slope value to test for sloping patterns

## Value

List containing two objects
