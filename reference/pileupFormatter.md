# Correctly formats pileup files.

Places columns in correct order and renames columns. Cleans the contig
labels to remove excess information after whitespace.

## Usage

``` r
pileupFormatter(pileup)
```

## Arguments

- pileup:

  A table containing contig names, read coverages averaged over 100 bp
  windows,and contig positions

## Value

dataframe
