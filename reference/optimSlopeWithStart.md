# Function for applying DIRECT global optimization algorithm

Set upper and lower limits for DIRECT and apply it to slope top, bottom,
and start values

## Usage

``` r
optimSlopeWithStart(
  minSlope,
  minSlopeSize,
  viralSubset,
  windowSize,
  leftOrRight,
  DirectMaxEval
)
```

## Arguments

- minSlope:

  The minimum slope value to test for sloping patterns

- minSlopeSize:

  The minimum width of sloping patterns.

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make.

## Value

List object
