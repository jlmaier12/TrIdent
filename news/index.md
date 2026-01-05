# Changelog

## Changes in 1.3.2

- Make the minimum sloping pattern widths a variable in
  TrIdentClassifier() (‘minSlopeSize’).
- Fixed minor bug where contigs smaller than a minimum pattern width
  would still allow one pattern variation to be tested.

## Changes in 1.3.1

- Make the median VLP:WC minimum ratio for HighCovNoPattern
  classifications a variable in TrIdentClassifier() (‘minHCNPRatio’).

## TrIdent V1.2 on Bioconductor 3.22!

## Changes in 1.1.5

- Added names to list items in pattern matching results.
- Fixed bug that caused an error when there were no “Prophage-like”,
  “Sloping”  
  and “HighCovNoPattern” classifications made.
- Added verbose for when there are no “Prophage-like”, “Sloping” and
  “HighCovNoPattern” classifications made.
- Cleaned specializedTransductionID() output table.

## Changes in 1.1.4

- Fixed bug caused by using the ‘specificContig’ parameter in the
  specializedTransductionID function.

## Changes in 1.1.3

- Added ability to normalize VLP-fraction:whole-community read coverage
  ratio by the number of reads in both the VLP-fraction and
  whole-community.

## Changes in 1.1.2

- Fixed bug where an error was caused when calculating a ratio with
  contigs in which the median whole-community non-prophage-like read
  coverage value was 0.
- Added ‘onlyPlot’ parameter in plotTrIdentResults() function.
- Added ‘logScale’ option in plotTrIdentResults() function.
- Removed ‘suggested filtering threshold’ feature.

## Changes in 1.1.1

- Fixed bug where an error was caused when calculating VLP:WC read
  coverage ratio with contigs in which the median whole-community read
  coverage value was

0.  

## Changes in 0.99.3

- Updated vignette to print a select few plots rather than all the
  output plots
- Added a text document to the inst/script/ directory with specific code
  and instructions for generating the sample pileup files used in the
  examples, vignette and README
- Fixed bug with verbose argument where verbose=FALSE still printed
  messages to the console
- Updated use of ‘if(verbose == TRUE)’ to ‘if(verbose)’

## Changes in 0.99.2

- Updated vignette and README with BiocManager installation instructions
- Moved Acknowledgements, Funding, and SessionInfo() to ‘supplemental
  information’ section in vignette
- Changed TOC_depth from 2 to 3 in vignette
- Updated formatting of NEWS file
- Grouped related functions in .R files for easier troubleshooting
- Add verbose = TRUE/FALSE argument to allow users to choose whether to
  print progress
- Added co-authors to DESCRIPTION file
- Added additional logic for input file validation
- Added better descriptions of input file requirements in the main
  function .man pages
- Added rendered summary histogram plot to output list
- Removed all usage of \<\<-
- Fixed improper usage of lapply

## Changes in 0.99.1

- Fixed unit tests for Bioc submission

## Changes in 0.99.0

- Bioc submission
