# gffSample

A subset of gene annotations for contigs from the whole-community. The
gff3 file was imported into R using the \`Bioc.gff:readGFF\` function.
Report...

## Usage

``` r
data('gffSample')
```

## Format

\## 'gffSample' A S4Vector::DataFrame with 1,085 rows and 30 columns:

- seqid:
- source:
- type:
- start:
- end:
- score:
- strand:
- phase:
- ID:
- genomedb_OC:
- genomedb_target:
- kegg_target:
- pfam_acc:
- pfam_desc:
- pfam_target:
- pfam_GO:
- tigrfam_acc:
- tigrfam_desc:
- tigrfam_target:
- eggnog_OC:
- eggnog_annot:
- eggnog_target:
- kegg_KO:
- kegg_desc:
- kegg_ecs:
- kegg_gos:
- kegg_pathwayNames:
- kegg_pathways:
- Name:
- rRNA_target:

## Source

\<https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5\>

## Details

This dataset represents the gene annotations obtained for the
whole-community contigs. The whole-community was generated from a
conventional mouse fecal homogenate.The whole-community DNA was
extracted and sequenced with Illumina (paired-end mode, 150 bp reads)
after which the metagenome was assembled. A subset of 10 contigs from
the pileup file were selected for this sample dataset. The contigs were
chosen because their associated read coverage patterns in the
VLP-fraction exemplify TrIdent's pattern-matching functionality across
classifications: NODE_617:Prophage-like, active/abundant, with spec
transduction NODE_135:Prophage-like, off one side of contig, no spec
transduction NODE_352:Sloping, left to right slope NODE_1088: Sloping,
right to left slope NODE_2060: Sloping, right to left slope with start
NODE_1401: None, no pattern match NODE_62: Prophage-like, with spec
transduction NODE_368: Prophage-like, not homogeneously
integrated/present, no spec transduction NODE_560: HighCovNoPattern
NODE_1165: None, filtered out To access the sequencing data used to
generate these gene annotations, refer to the reference below:
Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al.
Transductomics: sequencing-based detection and analysis of transduced
DNA in pure cultures and microbial communities. Microbiome 8, 158
(2020). https://doi.org/10.1186/s40168-020-00935-5
