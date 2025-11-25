# Specialized transduction search and plot

Search contigs classified as prophage-like for potential specialized
transduction and return the plot visualizing the search results.

## Usage

``` r
specTransductionSearch(
  contigName,
  VLPpileup,
  classifPatternMatches,
  classifSumm,
  windowSize,
  i,
  noReadCov,
  specTransLength,
  logScale
)
```

## Arguments

- contigName:

  The reference name of the contig currently being assessed (i.e
  "NODE_1")

- VLPpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping VLP-fraction
  reads to whole-community contigs

- classifPatternMatches:

  The pattern match information associated with each contig classified
  as prophage-like, sloping, or HighCovNoPattern

- classifSumm:

  The summary information associated with each contig classified as
  Prophage-like, Sloping, or HighCovNoPattern

- windowSize:

  The window size used to re-average read coverage pileups

- i:

  The index for the contig currently being assessed

- noReadCov:

  How many bp of no read coverage are encountered before searching
  stops? Default is 500.

- specTransLength:

  How many bp of read coverage to look for outside of prophage borders?
  Default is 2000.

- logScale:

  If TRUE, coverage is plotted in log10. If FALSE, raw coverage values
  are plotted. Default is FALSE.

## Value

List containing two objects
