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
TrIdent is a bioinformatics tool that automates the transductomics data analysis by automatically detecting, classifying and characterizing potential transducing events in a fraction of the time that it would take a manual-labeler. Transductomics is a DNA-sequencing based method for the detection and characterization of transduction events. Developed by Dr. Manuel Kleiner, transductomics relies on mapping reads from a virome (VLP-fraction) of a sample to contigs assembled from the metagenome (whole-community) of the same sample. Reads from bacterial DNA carried by viruses or VLPs (Viral-Like Particles) will map back to their bacterial contigs of origin creating read coverage patterns indicative of active transduction. 

Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al. Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities. Microbiome 8, 158 (2020). https://doi.org/10.1186/s40168-020-00935-5

VLP-fraction reads mapped to contigs from the whole-community may create the following read coverage patterns:
- Sloping pattern
  - Forms due to decreasing frequency of large DNA transfer moving away from a transduction initiation site.
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
    - When the ratio of read coverage on a contig between the VLP-fraction and the whole-community is high, its likely that the bacterial DNA is not in the VLP-fraction by random chance. I.e. contigs without sloping or block-like patterns but with high VLP-fraction:Whole-community read coverage ratios may represent the 'tails' to sloping patterns formed by generalized, lateral and gene transfer agent transducing events
    - Conversely, when the ratio of read coverage on a contig between the VLP-fraction and the whole-community is low, it's likely that any bacterial DNA in the VLP-fraction mapping to a contig is due to contamination and does not represent a transduction event.
 
With transductomics, a researcher can obtain information about the phage-host pairs involved in transduction, the types of transduction occuring and their locations for exploring associated gene annotations. Because transductomics detects DNA actively being carried by phage, it allows a glimpse into the ongoing transduction occuring in a microbiome rather than historical events. 

To obtain the data needed for transductomics, two sample types must be prepared for a microbiome of interest:
- Whole-community: Reprents the 'whole-community' (all bacteria, fungi, virus, etc) in the microbiome of interest
- VLP-fraction: Represents only the virus and 'viral-like particles' associated with the microbiome of interest
    - The VLP-fraction must be obtained via concentration of the viruses within your sample, then ultrapurification to remove bacterial cells and contaminating free bacterial DNA
 
After shotgun metagenomic sequencing of both samples (150bp short-read, paired-end, deep coverage), assemble a metagenome from the whole-community reads. Map the reads from both the VLP-fraction and whole-community to the metagenome contigs. Create pileup files using BBMap pileup.sh with 100bp windowsizes:
```{bash}
$ pileup.sh in=VLP_fraction_ReadMapping.bam out=VLP_Fraction.pileupcovstats bincov=VLP_fraction_pileup.bincov100 binsive=100 stdev=t
$ pileup.sh in=WholeCommunity_ReadMapping.bam out=WholeCommunity.pileupcovstats bincov=WholeCommunity_pileup.bincov100 binsive=100 stdev=t
```

Import the VLP_fraction_pileup.bincov100 and WholeCommunity_pileup.bincov100 files to R and run TrIdent to automatically detect, classify and characterize potential transduction events present on your contigs! 

## TrIdent description 

TrIdent first classifies contigs as 'Prophage-like', 'Gen/Lat/GTA', 'HighVLPWCRC', or 'None' using pattern-matching to detect patterns in read coverages.

- Prophage-like classifications include potential prophages or phage-inducible chromosomal islands (PICIs)
    - In addition to identifying prophages and PICIs on contigs, TrIdent also determines if they are highly active/abundant by assessing the associated read coverage in the whole-community fraction. If the read coverage of a prophage/PICI region is elevated above the non-prophage/PICI region of a contig, it is an indicator that the prophage/PICI is actively replicating or is highly abundant as its genome is represented in a higher ratio than its host bacteria.   
- Gen/Lat/GTA classifications include potential generalized, lateral or gene transfer agent (GTA) transduction events
- HighVLPWCRC classifications stand for 'High VLP-fraction:Whole-Community Read Coverage ratio' which means the contig has no pattern match but has an unusually high amount of bacterial DNA in the VLP-fraction and may represent the 'tail' of a sloping pattern formed by a Gen/Lat/GTA event.
- None classifications includes contigs with no pattern matches and low/no read coverage in the VLP-fraction.


Next, TrIdent can search contigs predicted as 'Prophage-like' for specialized transduction events to provide the user with a holistic view of the transduction actively occuring in their microbiome of interest. 

TrIdent is entirely reference-indepedent meaning that classifications made do not rely on the accuracy or availability of gene annotations. Prophage and PICIs that are not yet known or annotated can still be identified with TrIdent!

## Using TrIdent

**Import your pileup files:**
```{r}
VLP_fraction_pileup100 <- read.delim("Q:/VLP_fraction_pileup.bincov100", header=FALSE, comment.char="#")
WholeCommunity_pileup100 <- read.delim("Q:/WholeCommunity_pileup.bincov100", header=FALSE, comment.char="#")
```

**Run the TrIdent Classifier to classify contigs as Prophage-like, Gen/Lat/GTA, HighVLPWCRC, or None:** 
```{r}
Trident_results <- TrIdent_Classifier(VLP_fraction_pileup100, WholeCommunity_pileup100, windowsize=1000)

Summary_table <- Trident_results$Full_summary_table
Cleaned_table <- Trident_results$Cleaned_summary_table
Unused_contigs <- Trident_results$FilteredOut_contig_table
```

The output from TrIdent_Classifier is a list that contains five objects- 
1. Full_summary_table: A table containing the classifications and characterizations of all contigs that were not filtered out, **including** contigs that were classified as 'None'
2. Cleaned_summary_table: A table containing the classifications and characterizations of all contigs that were not filtered out, **excluding** contigs that were classified as 'None'
3. Pattern_MatchInfo: A list of information for each contig's pattern-match. This information is used by other functions in TrIdent. 
4. FilteredOut_contig_table: A table containing all contigs that were filtered out and the respective reason.
5. Windowsize: The windowsize used.


**Plot the results of the TrIdent_Classifier pattern-matching:**
```{r}
Trident_plots <- Plot_TrIdentPatternMatches(VLP_fraction_pileup100, WholeCommunity_pileup100, Trident_results)

#View all pattern matches to contigs classified as  Prophage-like, Gen/Lat/GTA, or HighVLPWCRC:
Trident_plots

#View/save specific plot:
Trident_plots$NODE_12
```

**Identify potential specialized transduction events on contigs classified as Prophage-like:**
```{r}
#Search all contigs classified as Prophage-like for specialized transduction
Spec_transduction <- SpecializedTransduction_ID(P_Spades3_100, Trident_results, noreadcov=500, spectranslength=2000)

Spec_transduction_summary <- Spec_transduction$Summary_table
Spec_transduction_plots <- Spec_transduction$Plots
Spec_transduction_Contig10 <- Spec_transduction_plots$NODE_10

#Search a specific contig classified as Prophage-like for specialized transduction
Spec_transduction_Contig1 <- SpecializedTransduction_ID(P_Spades3_100, Trident_results, specificcontig="NODE_1", noreadcov=500, spectranslength=2000)
```

SpecializedTransduction_ID can be used in two ways- it can search for specialized transduction in all contigs classified as Prophage-like or in a specific contig classified as Prophage-like. If you would like to search a specific contig for specialized transduction, include the contig name in quotes in the arguments (see example above). If you search a specific contig, the function will return the read coverage plot object for that contig. If you search all prophage-like contigs (aka you do not provide a specific contig's name), the function will return a list with the first object being a table summarizing the results of specialized transduction searching for each prophage-like contig. The second object is a list of read coverage plots for all prophage-like contigs. Each plot in the list is named by its reference name, i.e 'NODE_9'


If you would like to combine the output tables from TrIdent_Classifier and SpecializedTransduction_ID, try running this code:
```{r}
complete_trident_summary <- merge(Cleaned_table, Spec_transduction_summary, by="ref_name", all.x=TRUE)
```
