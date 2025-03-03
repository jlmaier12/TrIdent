## Changes in 0.99.4

+ Fixed bug where contigs in which the median whole-community read coverage 
value was 0 caused an error with VLP:WC read coverage ratio calculation

## Changes in 0.99.3

+ Updated vignette to print a select few plots rather than all the output 
plots
+ Added a text document to the inst/script/ directory with specific code and 
instructions for generating the sample pileup files used in the examples, vignette
and README
+ Fixed bug with verbose argument where verbose=FALSE still printed messages to 
the console
+ Updated use of 'if(verbose == TRUE)' to 'if(verbose)'

## Changes in 0.99.2

+ Updated vignette and README with BiocManager installation instructions
+ Moved Acknowledgements, Funding, and SessionInfo() to 'supplemental 
information'
section in vignette
+ Changed TOC_depth from 2 to 3 in vignette
+ Updated formatting of NEWS file
+ Grouped related functions in .R files for easier troubleshooting
+ Add verbose = TRUE/FALSE argument to allow users to choose whether to print 
progress
+ Added co-authors to DESCRIPTION file
+ Added additional logic for input file validation
+ Added better descriptions of input file requirements in the main function 
.man pages
+ Added rendered summary histogram plot to output list 
+ Removed all usage of <<-
+ Fixed improper usage of lapply

## Changes in 0.99.1

+ Fixed unit tests for Bioc submission

## Changes in 0.99.0

+ Bioc submission
