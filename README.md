# GenomicProcessingPipeline

Pipeline for processing "raw" genomic data (whole-genome sequencing, RNA sequencing, and target-capture sequencing)
_This code was used for data processing and analyses in Klabacka et al. (in prep)_

***

## Contents

-   [Project Background](#project-background)

***

## Project Backround

Studying factors driving natural variation in life-history strategies can help us understand how senescence evolves. Divergent ecotypes (slow-aging and fast-aging) of western terrestrial garter snake (Thamnophis elegans) provide a useful model system for examination of these factors. Here we examine gene expression and population genetics within and between these divergent ecotypes, and find support for hypothesized life-history divergence at the molecular level. We store our code for data processing and analyses, along with documentation for reproduction, within this repository.

***

## Study Design

#### Quantifying Gene Expression
32 individuals born and raised in the lab used within a 2 x 2 experimental design with heat treatment (27º C and 37º C) and ecotype (FA and SA) as variables.
#### Examining Targeted Sequence Variation
243 individuals genotyped for variant sites within 301 targeted genes (252 within networks of interest, 49 randomly selected)
Data for 94 of these individuals were sequenced using a target capture approach (Seq-Cap)
Data for 31 of these individuals were sequenced using a transcriptomic approach (RNA-Seq)
Data for 118 of these individuals were sequenced using a whole-genome approach (WGS)
We called SNPs using the Seq-Cap and RNA-Seq individuals, and then used this database to call SNPs from the same sites for the WGS individuals.

***

## Bioinformatics

### Bioinformatics Summary

Summarize bioinformatics here

### Gene Expression

Describe gene expression data processing and analyses here

### Sequence Variation

Scripts used for examination of targeted sequence variation are within the 'Examining-Sequence-Variation' directory. Here is a brief overview:

-   'q.main.sh' : This script executes functions from all other bash scripts for complete data processing.
-   'reads2vcf.sh' : This file contains functions for processing raw reads from RNA-Seq and Seq-Cap (cleaning, mapping, etc.) and calling SNPs 
-   'SNP_curation.sh' : This file contains functions for joining WGS data with Seq-Cap and RNA-Seq data, parsing SNPs into pairwise population files, inserting SNPs into multiple sequence alignments, calculating Tajima's D for each gene, and obtaining SNP proportions for their respective genes, transcripts, and coding regions.
-   'sift2vcf.sh' : This file contains functions for quantifying the functional implications of nonsynonymous SNPs and inserting these into a vcf.

