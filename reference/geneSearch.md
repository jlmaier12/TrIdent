# Search for gene annotations on positively classified contigs

Search contigs classified with TrIdent for gene-annotations that match a
provided key-word(s). Outputs read coverage plots for contigs with
matching annotations.

## Usage

``` r
geneSearch(
  TrIdentResults,
  VLPpileup,
  gff,
  searchCol,
  keyWords,
  inPatMat = FALSE,
  bpRange = 0,
  saveFilesTo,
  verbose = TRUE
)
```

## Arguments

- TrIdentResults:

  The output from \`TrIdent()\`.

- VLPpileup:

  A .txt file containing mapped sequencing read coverages averaged over
  100 bp windows/bins.

- gff:

  A .gff file imported using Bioc.gff::readGFF that is associated with
  the whole-community contigs

- searchCol:

  The exact column name (in quotes) in the gff file that you'd like to
  search for specific \`keyWords\`. Commonly annotation, product, gene,
  or acc columns, for example, however the specific name may change
  based on annotation tool used. For example: "pfam_desc".

- keyWords:

  The keyWord(s) to search for. Case independent. Searches will return
  the string that contains the matching keyWord. KeyWord(s) must be in
  quotes, comma-separated, and surrounded by c() i.e( c("antibiotic",
  "resistance", "drug") )

- inPatMat:

  TRUE or FALSE. If TRUE, only search for gene-annotations in the
  pattern-match region. Default is FALSE (i.e search the entire contig
  for the gene annotation key-words)

- bpRange:

  If \`inPatMat\` = TRUE, the user may specify the region (in base
  pairs) that should be searched to the left and right of the
  pattern-match region. Default is 0.

- saveFilesTo:

  Optional, Provide a path to the directory you wish to save output to.
  A folder will be made within the provided directory to store results.

- verbose:

  TRUE or FALSE. Print progress messages to console. Default is TRUE.

## Value

list of ggplot objects

## Examples

``` r
data("VLPFractionSamplePileup")
data("TrIdentSampleOutput")
data("gffSample")

patternMatches <- geneSearch(
  TrIdentResults = TrIdentSampleOutput,
  VLPpileup = VLPFractionSamplePileup,
  gff = gffSample,
  searchCol="pfam_desc",
  keyWords=c("toxin", "drug", "resistance", "phage")
)
#> Cleaning pileup file...
#> Searching for matching annotations...
#> 7 contigs have gene annotations that match one or more of the provided keyWords
```
