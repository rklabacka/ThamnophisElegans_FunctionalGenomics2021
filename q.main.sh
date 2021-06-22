#!/bin/sh

# This is the bash script used to process raw reads for the garter snake 
# targeted gene network study (for gene expression and populationi genetics)

# Prepared by Randy Klabacka

# -- Job Details for Hopper cluster at Auburn University-- #
#Give job a name
#PBS -N FullScript_May2020

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=1:ppn=10,walltime=05:00:00:00
#-- Indicate if\when you want to receive email about your job
#-- The directive below sends email if the job is (a) aborted, 
#-- when it (b) begins, and when it (e) ends
#PBS -m abe rlk0015@auburn.edu

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

source ./reads2vcf.sh
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
pythonScripts=/home/rlk0015/SeqCap/code/GenomicProcessingPipeline/pythonScripts

# prepare environment
loadModules
createWorkingEnvironment

## ------------------------
# begin raw read processing
## ------------------------
# Copy Sequence Capture raw reads into working environment
copyRawReadsDNA
# Quality check RNA reads
performFASTQC rawReadsDNA
## Note: RNA Seq reads were already copied
## and quality checked. If this needs to be
## done again, refer to the copyRawReadsDNA
## and performFastQC commands in reads2vcf.sh

# quality clean reads
performTrimmingPE DNA
performTrimmingPE RNA
performTrimmingSE RNA

# Quality check clean reads
performFASTQC cleanReadsDNA
performFASTQC cleanReadsRNA

# Copy reference genome into working environment
copyRef

# Map clean reads to reference genome
mapReadsDNA
mapPEReadsRNA
mapSEReadsRNA

# Change the sequence capture names to match sample ID
changeSeqCapNames

# Add sequencing read group information to each sample
readGroupsRNA 
readGroupsDNA

# Prep for SNP calling and calculate mapping stats
indexReference
prepForVariantCalling DNA
prepForVariantCalling RNA

### --------------------- Variant Calling ----------------------- ###
# Perform bqsr-'bootstraping', SNP calling, and hard filtration on Seq Cap data
cd $WorkingDirectory/GATKDNA
## Replicate 1
use-HaplotypeCaller 0 DNA
get-just-SNPs 0 
initial-VariantFiltration JustSNPs_0.vcf 0
use-BaseRecalibrator 0 DNA
use-BQSR 0 1 DNA
use-HaplotypeCaller 1 DNA
get-just-SNPs 1
use-BaseRecalibrator 1 DNA
use-AnalyzeCovariates 0 1 DNA
## -- Replicate 2
use-BQSR 1 2 DNA
use-HaplotypeCaller 2 DNA
get-just-SNPs 2
use-BaseRecalibrator 2 DNA
use-AnalyzeCovariates 1 2 DNA

# Perform bqsr-'bootstraping', SNP calling, and hard filtration on RNA-seq data
cd $WorkingDirectory/GATKRNA
## Replicate 1
use-HaplotypeCaller 0 RNA
get-just-SNPs 0
initial-VariantFiltration JustSNPs_0.vcf 0
use-BaseRecalibrator 0 RNA
use-BQSR 0 1 RNA
use-HaplotypeCaller 1 RNA
get-just-SNPs 1
use-BaseRecalibrator 1 RNA
use-AnalyzeCovariates 0 1 RNA
## -- Replicate 2
use-BQSR 1 2 RNA
use-HaplotypeCaller 2 RNA
get-just-SNPs 2
use-BaseRecalibrator 2 RNA
use-AnalyzeCovariates 1 2 RNA
## -- Merge RNA and DNA data
combine-VCF
cd $WorkingDirectory/variantFiltration 
## -- Eliminate potential RNA editing sites
removeRNAedits Merged
## -- Annotate variants
getNetworkFasta IILS
probes2gff Exons_2021.fa SeqCap
annotateVariants removedRNAedits SeqCap
## -- Initial Filter Variants
initial-VariantFiltration SeqCap_Annotated.vcf SeqCap_InitialFiltered
## -- Perform Hard Filtering
hard-VariantFiltration SeqCap_InitialFiltered SeqCap
## -- Get coding regions and exonic regions
getSpecificVariants SeqCap CDS
getSpecificVariants SeqCap Exons
## -- Copy final step of hard filtering (for nomenclature purposes)
cp SeqCap_HardFilterStep3.vcf SeqCap_Genes.vcf
## -- zip up everything
bgzip SeqCap_Genes.vcf
bgzip SeqCap_Exons.vcf
bgzip SeqCap_CDS.vcf
## -- plot the variants (you can do this with all of the files as desired)
gunzip Full_Exons.vcf.gz
plotVariants Full_Exons.vcf
bgzip Full_Exons.vcf

## -- files needed to be renamed, see reads2vcf.sh for details
renameSortedBAMs

# MAIN
loadModules
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
pythonScripts=/home/rlk0015/SeqCap/code/GenomicProcessingPipeline/pythonScripts
combineDatasets
sortVariants SeqCap_CDS   SeqCap_IndividualsToRemove 
sortVariants SeqCap_Exons  SeqCap_IndividualsToRemoves
sortVariants SeqCap_Genes  SeqCap_IndividualsToRemoves
sortVariants Full_CDS Full_IndividualsToRemove
sortVariants Full_Exons Full_IndividualsToRemove
sortVariants Full_Genes Full_IndividualsToRemove

getVariantBED SeqCap_CDS
getVariantBED SeqCap_Exons
getVariantBED SeqCap_Genes

and deleted... functionalAnnotation SeqCap
functionalAnnotation Full

getGeneVariants CDS
getGeneVariants Exons
getGeneVariants Genes
getGeneVariants CDS _missense
getGeneVariants CDS _synonymous

getTranscriptLengths CDS
getTranscriptLengths Exons
getTranscriptLengths Genes
getTranscriptLengths CDS _missense
getTranscriptLengths CDS _synonymous

vcf2faa
reference2faa
#+ COMPLETED moveCapturedGenes
#+ COMPLETED createMSA faa protein Sequences maskedMSA
#+ COMPLETED createMSA fna transcript Sequences maskedMSA

#+ COMPLETED createPopFiles
#+ COMPLETED createPairwiseVCFs
#+ COMPLETED getPairwisePopGen


