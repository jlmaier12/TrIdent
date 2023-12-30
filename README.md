![TrIDentLogo](https://github.com/jlmaier12/TrIdent/assets/45083046/15ef7ec7-49ac-48eb-86eb-d585d4b0869e)
**TrIdent-** **Tr**ansduction **Ident**ification: 
Automatic detection, classification and characterization of active transduction events in microbiomes. 

--------------------------------------------------------------------------------------------------------------------------------------------------------
## Thank you for being a beta-tester!
>**Some things to know:**
>
>- This page contains basic background information and basic usage of TrIdent. For more detailed information, visit the [tutorial here](./TrIdentTutorial.html).
>- TrIdent comes preloaded with a small sample dataset so users can follow along with the tutorial in their own R console. The tutorial takes ~5-10 minutes. 
>- Find TrIdent [installation instructions here](./TrIdentBeta_installation.html).
>- If you are using your own data, make sure it's correct! You need an ultrapurified virome and a metagenome from the same sample. Find detailed [input data requirements here](./TrIdentTutorial.html#input-data).
>- Please email me at jlmaier@ncsu.edu if you have any questions or issues.
>- If you try TrIdent to any extent, please let me know your thoughts [here](https://docs.google.com/forms/d/e/1FAIpQLSeGYRKpkhbBqWyioE0X_n6BoitlYjsy9SBM0GP6cOVOd3XOkA/viewform?usp=sf_link) so I can improve! 


## Background on Transductomics
TrIdent is a reference-independent bioinformatics tool that automates the transductomics data analysis by automatically detecting, classifying and characterizing potential transducing events. Transductomics is a DNA-sequencing based method for the detection and characterization of transduction events. Developed by Kleiner et al. (2020), transductomics relies on mapping reads from a virome (VLP-fraction) of a sample to contigs assembled from the metagenome (whole-community) of the same sample. Reads from bacterial DNA carried by viruses  and other VLPs (Virus-like particles) will map back to the bacterial contigs of origin creating read coverage patterns indicative of potential ongoing transduction.

Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al. Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities. Microbiome 8, 158 (2020). [https://doi.org/10.1186/s40168-020-00935-5]

VLP-fraction reads mapped to contigs from the whole-community may create the following read coverage patterns:
- Sloping pattern
  - Forms for transduction modes that lead to packaging of large regions of the host (bacterial) genome in a non-random manner.
  - May represent:
      - Generalized transduction
      - Lateral transduction
      - Gene transfer agents
- Block-like
    - Forms due to prophage or phage-inducible chromosomal island (PICI) reads mapping back to their integration site in a bacterial genome
    - Specialized transduction may appear as dense but low frequency read coverage on one or both sides of the block-like pattern. 
- No pattern
    - Forms when neither a block or sloping pattern is formed. 
    - Forms under several cirumstances:
        - No VLP-fraction reads map to the contig, i.e. no transduction
        - Contamination in the VLP-fraction causes a small amount of coverage across the contig
        - The 'tails' of sloping patterns formed by generalized, lateral and gene transfer agent transducing events
    - When the ratio of read coverage on a contig between the VLP-fraction and the whole-community is high, it is likely that the bacterial DNA is not in the VLP-fraction by random chance. I.e. contigs without sloping or block-like patterns but with high “VLP-fraction:Whole-community” read coverage ratios may represent the ‘tails’ to sloping patterns formed by generalized, lateral and gene transfer agent transducing events OR may represent unknown transduction pathways or contamination from other sources e.g. very small, very dense cells that were co-purified with VLPs.
    - Conversely, when the ratio of read coverage on a contig between the VLP-fraction and the whole-community is low, it is likely that any bacterial DNA in the VLP-fraction mapping to a contig is due to contamination and does not represent a transduction event.
 
With transductomics, a researcher can obtain information about the phage-host pairs involved (only for elements that integrate into host genome), the types of transduction occuring, and the region of the host genome that is potentially transduced, which allows exploration of transferred genes. To obtain the data needed for transductomics, a microbiome sample of interest is split to prepare two sub-sample types:
- Whole-community: Reprents the 'whole-community' (all bacteria, fungi, virus, etc) in the microbiome of interest
- VLP-fraction: Represents only the virus and 'viral-like particles' associated with the microbiome of interest
    - The VLP-fraction must be obtained by an appropriate ultra-purification protocol for your sample type to remove bacterial cells and contaminating free bacterial DNA.
 
After shotgun metagenomic sequencing of both samples (see Kleiner et al. (2020) for details on sequencing requirements), assemble a metagenome from the whole-community reads. Map the reads from both the VLP-fraction and whole-community to the metagenome contigs. Any read mapper can be used, however we recommend using [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/) with `ambiguous=random`, `qtrim=lr`, and `minid=0.97`. 

Create pileup files using `BBMap` pileup.sh with 100 bp window sizes:
```{bash}
$ pileup.sh in=YOUR_VLP_FRACTION_READ_MAPPING.bam out=VLP_Fraction.pileupcovstats bincov=VLP_fraction_pileup_bincov100.txt binsize=100 stdev=t
$ pileup.sh in=YOUR_WHOLE_COMMUNITY_READ_MAPPING.bam out=WholeCommunity.pileupcovstats bincov=WholeCommunity_pileup_bincov100.txt binsize=100 stdev=t
```

Import the VLP_fraction_pileup_bincov100.txt and WholeCommunity_pileup_bincov100.txt files to R and run TrIdent to automatically detect, classify and characterize potential transduction events present on your contigs! 

## TrIdent description 

TrIdent first classifies contigs as 'Prophage-like', 'Gen/Lat/GTA', 'HighVLPWCRC', or 'None' using pattern-matching to detect patterns in read coverages.

- Prophage-like classifications include potential prophages or phage-inducible chromosomal islands (PICIs)
    - In addition to identifying prophages and PICIs on contigs, TrIdent also determines if they are highly active/abundant by assessing the associated read coverage in the whole-community fraction. If the read coverage of a prophage/PICI region is elevated above the non-prophage/PICI region of a contig, it is an indicator that the prophage/PICI is actively replicating or is highly abundant as its genome is represented in a higher ratio than its host bacteria.  
- Gen/Lat/GTA classifications include potential generalized, lateral or gene transfer agent (GTA) transduction events
- HighVLPWCRC classifications stand for ‘High VLP-fraction:Whole-Community Read Coverage ratio’ which means the contig has no pattern match but has an unusually high amount of bacterial DNA in the VLP-fraction and may represent the ‘tail’ of a sloping pattern formed by a Gen/Lat/GTA event, unknown transduction pathways or contamination from other sources e.g. very small, very dense cells that were co-purified with VLPs.
- None classifications includes contigs with no pattern matches and low/no read coverage in the VLP-fraction.


Next, TrIdent can search contigs predicted as 'Prophage-like' for specialized transduction events to provide the user with a holistic view of the transduction actively occuring in their microbiome of interest. 

## Using TrIdent

**Import your pileup files:**
```{r}
VLP_fraction_pileup100 <- read.delim("VLP_fraction_pileup_bincov100.txt", header=FALSE, comment.char="#")
WholeCommunity_pileup100 <- read.delim("WholeCommunity_pileup_bincov100.txt", header=FALSE, comment.char="#")
```

**Run `TrIdent_Classifier()` to classify contigs as Prophage-like, Gen/Lat/GTA, HighVLPWCRC, or None:** 
```{r}
Trident_results <- TrIdent_Classifier(VLP_pileup=VLP_fraction_pileup100, WC_pileup=WholeCommunity_pileup100, minblocksize=10000, maxblocksize=200000, windowsize=1000)

Summary_table <- Trident_results$Full_summary_table
Cleaned_table <- Trident_results$Cleaned_summary_table
Unused_contigs <- Trident_results$FilteredOut_contig_table
```

The output from `TrIdent_Classifier()` is a list that contains five objects- 
1. Full_summary_table: A table containing the classifications and characterizations of all contigs that were not filtered out, **including** contigs that were classified as 'None'
2. Cleaned_summary_table: A table containing the classifications and characterizations of all contigs that were not filtered out, **excluding** contigs that were classified as 'None'
3. Pattern_MatchInfo: A list of information for each contig's pattern-match. This information is used by other functions in TrIdent. 
4. FilteredOut_contig_table: A table containing all contigs that were filtered out and the respective reason.
5. Windowsize: The windowsize used.

Note that `windowsize`, `minblocksize` and `maxblocksize` are user-defined variables. Please see the [TrIdent tutorial](./TrIdentTutorial.html) to understand how changing these variables may alter your results.

**Plot the results of the `TrIdent_Classifier()` pattern-matching:**
```{r}
Trident_plots <- Plot_TrIdentPatternMatches(VLP_pileup=VLP_fraction_pileup100, WC_pileup=WholeCommunity_pileup100, transductionclassifications=Trident_results)

#View all pattern matches to contigs classified as  Prophage-like, Gen/Lat/GTA, or HighVLPWCRC:
Trident_plots

#View/save specific plot:
Trident_plots$NODE_12
```

**Identify potential specialized transduction events on contigs classified as Prophage-like:**

`SpecializedTransduction_ID()` can be used in two ways- it can search for specialized transduction in all contigs classified as Prophage-like or in a specific contig classified as Prophage-like. If you would like to search a specific contig for specialized transduction, include the contig name in quotes in the arguments (see example above). If you search a specific contig, the function will return the read coverage plot object for that contig. If you search all prophage-like contigs (aka you do not provide a specific contig’s name), the function will return a list with the first object being a table summarizing the results of specialized transduction searching for each prophage-like contig. The second object is a list of read coverage plots for all prophage-like contigs. Each plot in the list is named by its reference name, e.g. ‘NODE_10’

Search all contigs classified as Prophage-like for specialized transduction
```{r}
Spec_transduction <- SpecializedTransduction_ID(VLP_pileup=P_Spades3_100, transductionclassifications=Trident_results, noreadcov=500, spectranslength=2000)

Spec_transduction_summary <- Spec_transduction$Summary_table
Spec_transduction_plots <- Spec_transduction$Plots
Spec_transduction_Contig10 <- Spec_transduction_plots$NODE_10
```

Search a specific contig classified as Prophage-like for specialized transduction
```{r}
Spec_transduction_Contig1 <- SpecializedTransduction_ID(VLP_pileup=P_Spades3_100, transductionclassifications=Trident_results, specificcontig="NODE_1", noreadcov=500, spectranslength=2000)
```
Note that `noreadcov` and `spectranslength` are user-defined variables. Please see the [TrIdent tutorial](./TrIdentTutorial.html) to understand how changing `noreadcov` and `spectranslength` may alter your results.

If you would like to combine the output tables from `TrIdent_Classifier()` and `SpecializedTransduction_ID()`, try this code:
```{r}
complete_trident_summary <- merge(Cleaned_table, Spec_transduction_summary, by="ref_name", all.x=TRUE)
```
