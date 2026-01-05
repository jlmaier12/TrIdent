# Sloping pattern translator

Translates a sloping pattern containing the initial jump-up in read
coverage across a contig. Translate the pattern 1000 bp at a time. Stop
translating when the pattern left on the contig reaches 20,000 bp.

## Usage

``` r
slopeTranslator(
  viralSubset,
  bestMatchInfo,
  windowSize,
  slopeChange,
  leftOrRight,
  minSlopeSize
)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- bestMatchInfo:

  The pattern-match information associated with the current best pattern
  match.

- windowSize:

  The window size used to re-average read coverage pileups

- slopeChange:

  A list containing pattern vector, slope value, and value of slope
  bottom

- leftOrRight:

  The direction of the sloping pattern. Either "Left" for left to right
  (neg) slopes or "Right" for right to left (pos) slopes.

- minSlopeSize:

  The minimum width of sloping patterns.

## Value

List
