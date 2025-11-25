# Specialized transduction plot

Plot search results of \`specializedTransductionID()\`

## Usage

``` r
specTransductionPlot(
  viralSubsetZoom,
  startPosBp,
  endPosBp,
  SpecTransLeft,
  specTransRight,
  contigName,
  classifPatternMatches,
  i,
  specTransSumm,
  logScale,
  classifSumm
)
```

## Arguments

- viralSubsetZoom:

  contig subset surrounding Prophage-like pattern-match

- startPosBp:

  Left border position

- endPosBp:

  Right border position

- SpecTransLeft:

  End position of spec transduction on left border

- specTransRight:

  End position of spec transduction on right border

- contigName:

  The reference name of the contig currently being assessed (i.e
  "NODE_1")

- classifPatternMatches:

  The pattern match information associated with each contig classified
  as prophage-like, sloping, or HighCovNoPattern

- i:

  The index for the contig currently being assessed

- specTransSumm:

  Results for spec transduction search

- logScale:

  If TRUE, coverage is plotted in log10. If FALSE, raw coverage values
  are plotted. Default is FALSE.

- classifSumm:

  The summary information associated with each contig classified as
  Prophage-like, Sloping, or HighCovNoPattern

## Value

ggplot object
