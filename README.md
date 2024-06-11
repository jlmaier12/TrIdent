![TrIDentLogo](https://github.com/jlmaier12/TrIdent/assets/45083046/15ef7ec7-49ac-48eb-86eb-d585d4b0869e)
**TrIdent-** **Tr**ansduction **Ident**ification: 
Automatic detection, classification and characterization of active transduction events in microbiomes. 

--------------------------------------------------------------------------------------------------------------------------------------------------------
## Thank you for being a beta-tester!
>**Some things to know:**
>
>- This page contains basic background information on TrIdent. To learn how TrIdent is used and the expected output, visit the [tutorial here](./TrIdentTutorial.html).
>- TrIdent comes preloaded with a small sample dataset so users can follow along with the tutorial in their own R console. The tutorial takes ~5-10 minutes. 
>- Find TrIdent [installation instructions here](./TrIdentBeta_installation.html).
>- If you are using your own data, make sure it's correct! You need an ultrapurified virome and a metagenome from the same sample. Find detailed [input data requirements here](./TrIdentTutorial.html#input-data).
>- Please email me at jlmaier@ncsu.edu if you have any questions or issues.
>- If you try TrIdent to any extent, please let me know your thoughts [here](https://docs.google.com/forms/d/e/1FAIpQLSeGYRKpkhbBqWyioE0X_n6BoitlYjsy9SBM0GP6cOVOd3XOkA/viewform?usp=sf_link) so I can improve! 
>- If read coverage pattern-matching is of interest to you, check out my other R package- [ProActive.](https://jlmaier12.github.io/ProActive/) 

## Background on Transductomics and TrIdent
TrIdent is a reference-independent bioinformatics tool that automates the transductomics data analysis by automatically detecting, classifying and characterizing potential transducing events. Transductomics is a DNA-sequencing based method for the detection and characterization of transduction events. Developed by Kleiner et al. (2020), transductomics relies on mapping reads from a virome (VLP-fraction) of a sample to contigs assembled from the metagenome (whole-community) of the same sample. Reads from bacterial DNA carried by viruses and other VLPs (Virus-like particles) will map back to the bacterial contigs of origin creating read coverage patterns indicative of potential ongoing transduction.

To obtain the data needed for transductomics, a microbiome sample of interest is split to prepare two sub-sample types:
- Whole-community: Represents the 'whole-community' (all bacteria, fungi, virus, etc) in the microbiome of interest
- VLP-fraction: Represents only the virus and 'viral-like particles' associated with the microbiome of interest
    - The VLP-fraction must be obtained by an appropriate ultra-purification protocol for your sample type to remove bacterial cells and contaminating free bacterial DNA.

Reference: Kleiner, M., Bushnell, B., Sanderson, K.E. et al. Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities. Microbiome 8, 158 (2020). [https://doi.org/10.1186/s40168-020-00935-5]
 
With transductomics and TrIdent, a researcher can obtain information about the phage-host pairs involved in transduction, the types of transduction occuring, and the region of the host genome that is potentially transduced, which allows exploration of transferred genes. 

