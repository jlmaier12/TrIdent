# Create histogram of normalized pattern-match scores

Plots a histogram of normalized match scores for all Prophage-like,
Sloping and HighCovNoPattern classifications and colors the plot based
on the classifications.

## Usage

``` r
resultsHisto(summaryList)
```

## Arguments

- summaryList:

  Classification summary table filtered to only include contigs with
  Prophage-like, Sloping and HighCovNoPattern classifications

## Value

ggplot object
