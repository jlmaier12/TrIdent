# Builds prophage-like block patterns

Build and translate a block pattern going off the left side, right side
and full length of the contig.

## Usage

``` r
blockBuilder(viralSubset, windowSize, minBlockSize, maxBlockSize)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

- minBlockSize:

  The minimum size of the prophage-like block pattern. Default is 10000
  bp.

- maxBlockSize:

  The maximum size of the prophage-like block pattern. Default is NA.

## Value

List containing three objects
