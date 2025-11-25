# Prophage-like border finder

Find borders of Prophage-like patterns with more specificity than
pattern-matching using 100 bp window pileups and sliding standard
deviation technique.

## Usage

``` r
prophageLikeBorders(viralSubset, classificationPatterns, i, windowSize)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- classificationPatterns:

  The pattern match information associated with each contig classified
  as Prophage-like, Sloping, or HighCovNoPattern

- i:

  The index for the contig currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

## Value

List
