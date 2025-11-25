# Change slope of sloping pattern

Change the value of the slope used for the sloping pattern-match

## Usage

``` r
changeSlope(
  leftOrRight,
  slopeBottom,
  halfToMaxReadCov,
  cov,
  viralSubset,
  windowSize
)
```

## Arguments

- leftOrRight:

  Generate pattern for negative slope (left to right, i.e. 'Left') or
  positive slope (right to left, i.e. 'Right')

- slopeBottom:

  The value for the bottom of the sloping value

- halfToMaxReadCov:

  Half of the max VLP-fraction read coverage divided by 10

- cov:

  The value for the top of the slope

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileup

## Value

List
