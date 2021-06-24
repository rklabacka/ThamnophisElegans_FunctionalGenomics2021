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

#### Scripts & Files
Scripts and coding files used for examination of targeted sequence variation are within the 'Examining-Sequence-Variation' directory. Here is a brief overview:

-   'q.main.sh' : This script executes functions from all other bash scripts for complete data processing.
-   'reads2vcf.sh' : This file contains functions for processing raw reads from RNA-Seq and Seq-Cap (cleaning, mapping, etc.) and calling SNPs 
-   'SNP_curation.sh' : This file contains functions for joining WGS data with Seq-Cap and RNA-Seq data, parsing SNPs into pairwise population files, inserting SNPs into multiple sequence alignments, calculating Tajima's D for each gene, and obtaining SNP proportions for their respective genes, transcripts, and coding regions.
-   'sift2vcf.sh' : This file contains functions for quantifying the functional implications of nonsynonymous SNPs and inserting these into a vcf.

#### Workflow
Bioinformatic pipelines can be complex and complicated. Here I will describe the general workflow, providing descriptions where some detail is necessary. For a more-detailed description, reading through the scripts/files themselves (and potentially documentation for the tools/packages used) may be necessary.

1.  Raw reads to mapped alignment

    We begin with raw '.fastq' files which we received from the genomic sequencing company. We need to clean these reads to (A) remove the adapter sequence and (B) remove low-quality information that may be incorrect due to sequencing error. To do this, we first check the quality using the program [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). This program provides information about our reads, including position-specific quality scores, read-wide quality scores, and adapter content. Here is an example of the position quality scores for our Seq-Cap reads: 
![Raw Read FastQC Quality](./Examining-Sequence-Variation/images/RawReadsFastQC.png)

    You'll notice that the quality of each base call ([Phred score](https://en.wikipedia.org/wiki/Phred_quality_score)] decreases toward the end of the reads. The FastQC output (an html file) contains many other plots to help assess read quality. To clean up our reads and remove the sequence adapters (which were used for binding, sequence initiation, and indexing), we use the program [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). This program will clean up our reads using our specified parameters.

    Following read cleaning, we check the read quality again. This time our position quality scores for our Seq-Cap reads look much better:
![Clean Read FastQC Quality](./Examining-Sequence-Variation/images/CleanReadsFastQC.png)

    After cleaning our reads, we are ready to map them to a reference. This can be challenging from a study design perspective; the decision for how to map can be a tricky one. If you have a reference genome for your focal taxon (which we luckily did), you can simply map to this. Alternatively, you can map to a transcriptome or the genome of a closely-related species. We map our cleaned reads using two approaches: (1) for our reads from Seq-Cap, we mapped using the program [BWA](https://hpc.nih.gov/apps/bwa.html), (2) for our reads from RNA-Seq, we used [HiSat2](http://daehwankimlab.github.io/hisat2/). The approach you take depends on your nucleotide type and sequencing approach (e.g., reads from single-end sequencing should be treated differently than those from paired-end sequencing). Mapping will use the clean .fastq files to create a [.sam (sequence alignment map)](https://en.wikipedia.org/wiki/SAM_(file_format)) file. These can be converted to the compressed version, .bam, using [Samtools](https://www.htslib.org/) to increase downstream processing efficiency.

2.  Mapped Alignment to Variant Calls

    Once reads have been mapped and stored in an alignment file, variation at specific positions can be identified. To prepare for variant calling, it is important to identify reads that might be duplicates (e.g., from PCR). We mark these duplicates using the MarkDuplicates tool from [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-).

    After marking duplicates, we perform a round of variant calling. For this project, we are specifically interested in single nucleotide polymorphisms (SNPs) that we can identify using various tools within the software package [GATK](https://gatk.broadinstitute.org). We follow the [GATK best practices workflow for germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) and implement suggestions regarding [base score recalibration](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) in non-model organisms. Like many software packages, GATK is updated regularly. We used version [4.1.7](https://gatk.broadinstitute.org/hc/en-us/articles/360042912311--Tool-Documentation-Index), but advise others to be aware that changes in versions may affect functionality/outcomes.

    First, we call variants using [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller). We then create a variant database and select only SNPs using [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360042477052-GenomicsDBImport), [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360042914991-GenotypeGVCFs), and [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360042913071-SelectVariants). We then use [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360042477652-VariantFiltration) to select SNPs we identify with high confidence.

    Next, we perform what GATK developers refer to as "bootstrapping." We use our high-confidence SNPs to generate a recalibration table (table_0) with our .bam file using [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042477672-BaseRecalibrator). We then perform base score recalibration to create a new .bam file using [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360042476852-ApplyBQSR). With this new .bam file, we repeat the variant calling process described above (with the exception of variant filtration) and create another recalibration table (table_1). We then compare table_0 with table_1 using [AnalyzeCovariates](https://gatk.broadinstitute.org/hc/en-us/articles/360042911971-AnalyzeCovariates). This tool outputs a pdf that compares assigned quality scores and their accuracy between table_0 and table_1 within each individual. Here is an example of the exported plots:
![Clean Read FastQC Quality](./Examining-Sequence-Variation/images/AnalyzeCovariates_0.png)
