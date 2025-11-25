# Full block pattern-translator

Translates full block-pattern across a contig. Translate the pattern
1000 bp at a time. Stop translating when the pattern is 5000 bp from the
end of the contig.

## Usage

``` r
blockTranslator(viralSubset, bestMatchInfo, windowSize, pattern)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- bestMatchInfo:

  The information associated with the current best pattern-match.

- windowSize:

  The window size used to re-average read coverage pileups

- pattern:

  A vector containing the values associated with the block pattern

## Value

List
