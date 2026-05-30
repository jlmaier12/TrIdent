# DIRECT block pattern optimizer

Optimizes prophage-like block patterns using DIRECT global optimization.

## Usage

``` r
directBlockBuilder(
  viralSubset,
  windowSize,
  minBlockSize,
  blockLength,
  blockLengthFull,
  minReadCov,
  maxReadCov,
  startingCovs,
  DirectMaxEval
)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

- minBlockSize:

  The minimum size of the prophage-like block pattern.

- blockLength:

  Maximum left and right block pattern length

- blockLengthFull:

  Maximum full block pattern length

- minReadCov:

  Baseline value outside of the block

- maxReadCov:

  Maximum VLP-fraction read coverage value

- startingCovs:

  Candidate coverage values used to initialize the search

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make. Default is 100.

## Value

List containing three objects
