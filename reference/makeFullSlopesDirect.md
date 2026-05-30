# Make full slope patterns

Makes slope patterns sloping either left to right (Left) or right to
left (right) across the contig being assessed using DIRECT dimensions

## Usage

``` r
makeFullSlopesDirect(leftOrRight, viralSubset, dims, windowSize, minSlope)
```

## Arguments

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- dims:

  Slope top and bottom values

- windowSize:

  The window size used to re-average read coverage pileups

- minSlope:

  The minimum slope value to test for sloping patterns

## Value

List
