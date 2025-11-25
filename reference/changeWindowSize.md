# Change the read coverage rolling mean window size

Re-averages window sizes of read coverage averages. Start with 100bp
windows always. Cannot make window size less than 100bp.

## Usage

``` r
changeWindowSize(cleanPileup, windowSize)
```

## Arguments

- cleanPileup:

  A read coverage dataset that has been cleaned and reformatted.

- windowSize:

  The number of base pairs to average coverage values over. Options are
  100, 500, 1000, or 2000 only!

## Value

Dataframe
