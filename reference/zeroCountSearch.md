# Counts zero values to the left and right of prophage-like borders

Checks to see at which point the number of consecutive zero values to
the left and right of the prophage-like pattern match borders equals the
noReadCov parameter

## Usage

``` r
zeroCountSearch(startOrEnd, viralSubsetZoom, startOrEndPosRow, noReadCov)
```

## Arguments

- startOrEnd:

  searching the start (left side) or end (right side) of the
  prophage-like pattern-match

- viralSubsetZoom:

  viralSubset dataframe subsetted to 50,000 bp outside the pattern match
  borders

- startOrEndPosRow:

  The row index of the start or end position of the prophage-like
  pattern match

- noReadCov:

  How many bp of no read coverage are encountered before specialized
  transduction searching stops? Default is 500.

## Value

List
