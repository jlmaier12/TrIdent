# No pattern pattern-match

A horizontal line at the mean or median coverage should be an optimal
pattern-match if the contig read coverage displays no sloping or block
patterns

## Usage

``` r
noPattern(viralSubset, searchMethod, DirectMaxEval, globalLocal)
```

## Arguments

- viralSubset:

  A subset of the read coverage pileup that pertains only to the contig
  currently being assessed

- searchMethod:

  Search method to use. Either "grid" for the original grid search or
  "direct" for DIRECT global optimization.

- DirectMaxEval:

  Maximum number of DIRECT evaluations to make.

- globalLocal:

  Use global or local DIRECT search. Default is local.

## Value

List
