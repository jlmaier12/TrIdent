# TrIdentSampleOutput

The TrIdentClassifier output from the VLPFractionSamplePileup and
WholeCommunitySamplePileup files run with default parameters Report...

## Usage

``` r
data('TrIdentSampleOutput')
```

## Format

\## 'TrIdentSampleOutput' A list with 6 objects:

- SummaryTable:

  A dataframe containing classifications for all contigs that were
  processed with pattern-matching

- CleanedSummaryTable:

  SummaryTable dataframe filtered to remove contigs that recieved a
  'None' classification

- PatternMatchInfo:

  A list of lists containing pattern-match information for each
  classified contig

- FilteredOutContigTable:

  A dataframe containing names of contigs that were filtered out prior
  to pattern-matching

- windowSize:

  windowSize used in TrIdentClassifier function (1000)

- ResultHistogram:

  a histogram displaying the overall abundance and quality of
  pattern-matches in addition to the composition of classifications. The
  displayed pattern-match scores are normalized by dividing each score
  by its associated contig length. The scores are normalized to
  visualize the overall quality of pattern-matching for the entire
  dataset.

## Details

A list object produced by the TrIdentClassifier function run on the
VLPFractionSamplePileup and WholeCommunitySamplePileup files run with
default parameters
