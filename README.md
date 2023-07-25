# TrIdent
**TrIdent-** **Tr**ansduction **Ident**ification: 
Automatic detection, classification and characterization of active transduction events in microbiomes. 

## Background on Transductomics
TrIdent is a bioinformatics tool that automates the transductomics data analyis by automatically detecting, classifying and characterizing potential transducing events in a fraction of the time that it would take a manual-labeler. Transductomics is a DNA-sequencing based method for the detection and characterization of transduction events. Developed by Dr. Manuel Kleiner, the method relies on the mapping reads from a virome(VLP_fracion) of a sample to contigs assembled from the metagenome(whole-community) of the same sample. Reads from bacterial DNA carried by viruses or VLPs (Viral-Like Particles) will map back to their bacterial contigs of origin creating read coverage patterns indicative of active transduction. 

VLP-fraction reads mapped to contigs from the whole-community will create the following read coverage patterns:
- Sloping pattern
  - Forms due to decreasing frequency of large DNA transfer startin from a transduction initiation site.
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
 
After shotgun metagenomic sequencing of both samples, assemble a metagenome from the whole-community reads. Map the reads from both the VLP-fraction and whole-community to the metagenome contigs. Create pileup files using BBMap pileup.sh with 100bp windowsizes:
`{bash}
$ pileup.sh in=VLP_fraction_ReadMapping.bam out=VLP_Fraction.pileupcovstats bincov=VLP_fraction_pileup.bincov100 binsive=100 stdev=t
$ pileup.sh in=WholeCommunity_ReadMapping.bam out=WholeCommunity.pileupcovstats bincov=WholeCommunity_pileup.bincov100 binsive=100 stdev=t
`

Import the VLP_fraction_pileup.bincov100 and WholeCommunity_pileup.bincov100 files to R and run TrIdent to automatically detect, classify and characterize potential transduction events present on your contigs! 

## TrIdent description 

TrIdent first classifies contigs as 'Prophage-like', 'Gen/Lat/GTA', 'HighVLPWCRC', or 'None' using pattern-matching to detect patterns in read coverages.

- Prophage-like classifications include potential prophages or phage-inducible chromosomal islands (PICIs)
    - In addition to identifying prophages and PICIs on contigs, TrIdent also determines if they are highly active/abundant by assessing the associated read coverage in the whole-community fraction. If the read coverage of a prophage/PICI region is elevated above the non-prophage/PICI region of a contig, it is an indicator that the prophage/PICI is actively replicating or is highly abundant as its genome is represented in a higher frequency than its host bacteria.   
- Gen/Lat/GTA classifications include potential generalized, lateral or gene transfer agent (GTA) transduction events
- HighVLPWCRC classifications stand for 'High VLP-fraction:Whole-Community Read Coverage ratio' which means the contig does not have a pattern match but has an unusually high amount of bacterial DNA in the VLP-fraction and may represent the 'tail' of a sloping pattern formed by a Gen/Lat/GTA event.
- None classifications includes contigs with no pattern matches and low/no read coverage in the VLP-fraction.


Next, TrIdent can search contigs predicted as 'Prophage-like' for potential specialized transduction events to provide the user with a holistic view of the transduction actively occuring in their microbiome of interest. 

TrIdent is entirely reference-indepedent meaning that classifications made do not rely on the accuracy or availability of gene annotations. Prophage and PICIs that are not yet known or annotated can still be identified with TrIdent!

## Using TrIdent

Run the TrIdent Classifier to classify contigs as prophage-like
`{R}

`


