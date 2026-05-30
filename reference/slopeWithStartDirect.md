# Sloping pattern with an initial jump-up in read coverage

Build, translate, and change slope of sloping pattern with slope start
using DIRECT dimensions

## Usage

``` r
slopeWithStartDirect(
  viralSubset,
  windowSize,
  minSlope,
  minSlopeSize,
  leftOrRight,
  dims
)
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

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- dims:

  Slope top, bottom, and start position

## Value

List object
