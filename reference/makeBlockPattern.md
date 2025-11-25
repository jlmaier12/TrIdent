# Make block patterns for pattern-matching

Make full, left and right block patterns for Prophage-like
classifications

## Usage

``` r
makeBlockPattern(
  viralSubset,
  windowSize,
  fullLeftRight,
  blockLength,
  nonBlock,
  minReadCov,
  cov
)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- windowSize:

  The window size used to re-average read coverage pileups

- fullLeftRight:

  The block pattern variation being built

- blockLength:

  Maximum block pattern length

- nonBlock:

  Maximum non-block pattern length

- minReadCov:

  Either 0 or the minimum VLP-fraction read coverage value

- cov:

  The height value of the block pattern

## Value

List containing two objects
