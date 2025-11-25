# Pattern-builder

Builds the pattern (vector) associated with the best pattern-match' for
each contig classified as Prophage-like, Sloping, or HighCovNoPattern.

## Usage

``` r
patternBuilder(viralSubset, classifList, classification, rowIndex)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- classifList:

  A list containing pattern match information associated with all
  classified contigs.

- classification:

  The contig's classification assigned by the TrIdentClassifier function

- rowIndex:

  The list index associated with each contig's pattern-match information

## Value

Vector
