# Pattern-match size calculator

Calculate the size (bp) of the matching region for Prophage-like and
Sloping patterns

## Usage

``` r
patternMatchSize(classifSumm, classifList, windowSize, verbose)
```

## Arguments

- classifSumm:

  Classification summary table

- classifList:

  A list containing pattern match information associated with all contig
  classifications

- windowSize:

  The window size used to re-average read coverage pileups

- verbose:

  TRUE or FALSE. Print progress messages to console. Default is TRUE.

## Value

dataframe
