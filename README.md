![TrIDentLogo](https://github.com/jlmaier12/TrIdent/assets/45083046/15ef7ec7-49ac-48eb-86eb-d585d4b0869e)
**TrIdent-** **Tr**ansduction **Ident**ification: 
Automatic detection, classification and characterization of active transduction events in microbiomes. 

--------------------------------------------------------------------------------------------------------------------------------------------------------
## Install:
```{r}
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("jlmaier12/TrIdent")
library(TrIdent)
```
## Background on Transductomics and TrIdent
TrIdent automates the analysis of transductomics data by detecting, classifying, and characterizing read coverage patterns associated with potential transduction events. Transductomics, developed by Kleiner et al. (2020), is a DNA sequencing-based method for the detection and characterization of transduction events in pure cultures and complex communities. Transductomics relies on mapping sequencing reads from a viral-like particle (VLP)-fraction of a sample to contigs assembled from the metagenome (whole-community) of the same sample. Reads from bacterial DNA carried by VLPs will map back to the bacterial contigs of origin creating read coverage patterns indicative of ongoing transduction. **The read coverage patterns detected represent DNA being actively carried or transduced by VLPs. The read coverage patterns do not represent complete transduction events (i.e integration of transduced DNA into new bacterial chromosomes).** 

To obtain the data needed for transductomics, a microbiome sample of interest is split to prepare two sub-sample types:
- Whole-community: Represents the 'whole-community' (all bacteria, fungi, virus, etc) in the microbiome of interest
- VLP-fraction: Represents only the virus and 'viral-like particles' associated with the microbiome of interest
    - The VLP-fraction must be obtained by an appropriate ultra-purification protocol for your sample type to remove bacterial cells and contaminating free bacterial DNA.

With transductomics and TrIdent, a researcher can obtain information about the phage-host pairs involved in transduction, the types of transduction occuring, and the region of the host genome that is potentially transduced, which allows exploration of transferred genes. 

**Reference:** Kleiner, M., Bushnell, B., Sanderson, K.E. et al. Transductomics: sequencing-based detection and analysis of transduced DNA in pure cultures and microbial communities. Microbiome 8, 158 (2020). <https://doi.org/10.1186/s40168-020-00935-5>

## Quick-start
```{r}
## Run first
TrIdentOutput <- TrIdentClassifier(VLPpileup=VLPFractionSamplePileup, 
                                   WCpileup=WholeCommunitySamplePileup)

## Run second
TrIdentPlots <- plotTrIdentResults(VLPpileup=VLPFractionSamplePileup, 
                                   WCpileup=WholeCommunitySamplePileup, 
                                   TrIdentResults=TrIdentOutput)

## Run third
SpecTransduc <- specializedTransductionID(VLPpileup=VLPFractionSamplePileup, 
                                          TrIdentResults=TrIdentOutput)
```

### Input files

TrIdent detects read coverage patterns using a pattern-matching algorithm that operates on pileup files. A pileup file is a file format where each row summarizes the 'pileup' of reads at specific genomic locations. Pileup files can be used to generate a rolling mean of read coverages and associated base pair positions across a metagenome assembly which reduces data size while preserving read coverage patterns. **TrIdent requires that input pileups files be generated using a 100 bp window/bin size.**  

Some read mappers, like [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/), will allow for the generation of pileup files in the [`bbmap.sh`](https://github.com/BioInfoTools/BBMap/blob/master/sh/bbmap.sh) command with the use of the `bincov` output with the `covbinsize=100` parameter/argument. **Otherwise, BBMap's [`pileup.sh`](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh) can convert .bam files produced by any read mapper to pileup files compatible with TrIdent using the `bincov` output with `binsize=100`.**

TrIdent requires two pileup files from a transductomics dataset as input:

* A VLP-fraction pileup: Sequencing reads from a sample's ultra-purified VLP-fraction mapped to the whole-community metagenome assembly from the same sample.
* A whole-community pileup: Sequencing reads from a sample's whole-community mapped to the whole-community metagenome from the same sample. 

**The data used for each pileup file must originate from the same sample. Pileup files must use a 100 bp window/bin size for the rolling mean.** 

Transductomics sample preparation, sequencing procedures, and analysis methods are detailed in [Kleiner et al. (2020)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00935-5)

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

If read coverage pattern-matching is of interest to you, check out my other R package- [ProActive.](https://jlmaier12.github.io/ProActive/) 
