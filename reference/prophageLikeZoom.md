# Prophage-like pattern zoom

'Zoom-in' on (aka subset) desired region surrounding block pattern.

## Usage

``` r
prophageLikeZoom(viralSubset, classificationPatterns, i, zoom, windowSize)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- classificationPatterns:

  The pattern match information associated with each contig classified
  as Prophage-like, sloping, or HighCovNoPattern

- i:

  The index for the contig currently being assessed

- zoom:

  The number of rows outside the start and stop positions of the block
  pattern to zoom-in on

- windowSize:

  The window size used to re-average read coverage pileups

## Value

Dataframe
