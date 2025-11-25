# Determine Prophage-like read coverage elevation in whole-community

Determines whether a detected Prophage-like genetic element has read
coverage in the whole-community that is either elevated or depressed
compared to the average read coverage of the non-prophage region.

## Usage

``` r
prophageLikeElevation(
  classifSummTable,
  prophageLikeClassifList,
  VLPpileup,
  WCpileup,
  windowSize,
  verbose
)
```

## Arguments

- classifSummTable:

  Classification summary table

- prophageLikeClassifList:

  A list containing pattern match information associated with all
  contigs classified as Prophage-like.

- VLPpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping VLP-fraction
  reads to whole-community contigs

- WCpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping whole-community
  reads to whole-community contigs

- windowSize:

  The window size used to re-average read coverage pileups

- verbose:

  TRUE or FALSE. Print progress messages to console. Default is TRUE.

## Value

dataframe
