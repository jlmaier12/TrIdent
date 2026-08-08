# Function for applying DIRECT global optimization algorithm

Set upper and lower limits for DIRECT and apply it to the noPattern
pattern

## Usage

``` r
optimNoPattern(viralSubset, DirectMaxEval, globalLocal)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make.

- globalLocal:

  Use global or local DIRECT search. Default is local.

## Value

List
