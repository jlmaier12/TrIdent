# Sloping pattern with an initial jump-up in read coverage

Build, translate, and change slope of sloping pattern with slope start

## Usage

``` r
slopeWithStart(viralSubset, windowSize, minSlope, minSlopeSize)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

- minSlope:

  The minimum slope value to test for sloping patterns

- minSlopeSize:

  The minimum width of sloping patterns.

## Value

List containing two objects
