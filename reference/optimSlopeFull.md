# Function for applying DIRECT global optimization algorithm

Set upper and lower limits for DIRECT and apply it to slope top and
bottom

## Usage

``` r
optimSlopeFull(
  minSlope,
  viralSubset,
  windowSize,
  leftOrRight,
  DirectMaxEval,
  globalLocal
)
```

## Arguments

- minSlope:

  The minimum slope value to test for sloping patterns

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileup

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make.

- globalLocal:

  Use global or local DIRECT search. Default is local.

## Value

List
