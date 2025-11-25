# VLP-fraction:whole-community read coverage ratio calculator

Calculate the VLP-fraction:whole-community read coverage ratio for every
contig using the median read coverage values. If the ratio is greater
than 2 (i.e VLP-fraction read coverage is, on average, at least double
the whole-community read coverage), then the contig is classified as
HighCovNoPattern. If the number of VLP-fraction and whole-community
reads used for mapping are provided, then the VLP/WC ratio value will be
normalized to the sizes of the VLP and WC read sets.

## Usage

``` r
VLPtoWCRatioCalc(
  classifSumm,
  WCpileup,
  VLPpileup,
  VLPReads,
  WCReads,
  minHCNPRatio
)
```

## Arguments

- classifSumm:

  Classification summary table

- WCpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping whole-community
  reads to whole-community contigs

- VLPpileup:

  A table containing contig names, coverages averaged over 100 bp
  windows, and contig positions associated with mapping VLP-fraction
  reads to whole-community contigs

- VLPReads:

  The number of VLP-fraction reads used for mapping and creation of
  pileup.

- WCReads:

  The number of WC reads used for mapping and creation of pileup.

- minHCNPRatio:

  The minimum VLP:WC ratio value used for HighCovNoPattern
  classifications. Default is 2. (i.e the median VLP-fraction coverage
  must be at least 2x the median WC read coverage to be classified as
  HighCovNoPattern).

## Value

dataframe
