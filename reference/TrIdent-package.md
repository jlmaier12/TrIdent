# **TrIdent** - **Tr**ansduction **Ident**ification

Automatic detection, classification and characterization of transduction
events in transductomics datasets using read coverage pattern-matching.

Please see \[Transductomics: sequencing-based detection and analysis of
transduced DNA in pure cultures and microbial communities\] (
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5)
for more information on the transductomics method, data and analysis
workflow.

## Details

The three main functions in TrIdent are:

1.  [`TrIdentClassifier`](https://jlmaier12.github.io/TrIdent/reference/TrIdentClassifier.md)
    performs the pattern-matching, classification and characterization
    of read coverage patterns on contigs.

2.  [`plotTrIdentResults`](https://jlmaier12.github.io/TrIdent/reference/plotTrIdentResults.md)
    plots the results from
    [`TrIdentClassifier()`](https://jlmaier12.github.io/TrIdent/reference/TrIdentClassifier.md)

3.  [`specializedTransductionID`](https://jlmaier12.github.io/TrIdent/reference/specializedTransductionID.md)
    searches contigs classified as Prophage-like by
    [`TrIdentClassifier()`](https://jlmaier12.github.io/TrIdent/reference/TrIdentClassifier.md)
    for potential specialized transduction

## See also

Useful links:

- <https://github.com/jlmaier12/TrIdent>

- <https://jlmaier12.github.io/TrIdent/>

- Report bugs at <https://github.com/jlmaier12/TrIdent/issues>

## Author

Jessie Maier <jlmaier@ncsu.edu> & Jorden Rabasco <jrabasc@ncsu.edu>
