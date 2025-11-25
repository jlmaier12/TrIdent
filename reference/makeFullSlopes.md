# Make full slope patterns

Makes slope patterns sloping either left to right (Left) or right to
left (right) across the contig being assessed.

## Usage

``` r
makeFullSlopes(leftOrRight, viralSubset, newMax, minReadCov, windowSize)
```

## Arguments

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- newMax:

  A value for the top of the sloping pattern that is slightly higher
  than the maximum coverage value on the viralSubset

- minReadCov:

  Minimum read coverage value of the viralSubset

- windowSize:

  The window size used to re-average read coverage pileups

## Value

List
