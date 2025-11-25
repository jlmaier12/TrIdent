# Translate left and right block patterns across contig

Translates left and right block patterns across contigs 1000 bp at a
time

## Usage

``` r
leftRightBlockTranslater(
  viralSubset,
  pattern,
  leftOrRight,
  windowSize,
  minReadCov,
  cov,
  bestMatchInfo,
  minBlockSize
)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- pattern:

  The pattern vector being translated

- leftOrRight:

  Is the left or right block pattern being translated

- windowSize:

  The window size used to re-average read coverage pileups

- minReadCov:

  The baseline value used for the region outside of the block pattern
  (either 0 or the minimum VLP-fraction read coverage for the contig)

- cov:

  The height value currently being used for the block pattern

- bestMatchInfo:

  The information associated with the current best pattern-match.

- minBlockSize:

  The minimum size of the Prophage-like block pattern. Default is 10,000
  bp.

## Value

List
