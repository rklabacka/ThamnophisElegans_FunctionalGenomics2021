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

WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
pythonScripts=/home/rlk0015/SeqCap/code/GenomicProcessingPipeline/pythonScripts

function loadModules {
  module load fastqc/11.5
  module load gnu-parallel/20160322
  module load trimmomatic/0.36
  module load bwa/0.7.15
  module load samtools/1.3.1
  module load picard/2.4.1
  module load hybpiper/1
  module load xz/5.2.2
  module load htslib/1.3.3
  module load python/3.6.4
  module load gatk/4.1.7.0
  module load R/3.6.0
  module load bedtools/2.29.0
  module load ncbi-blast/2.9.0+
  module load bcftools/1.3.2
  module load vcftools/v0.1.17
  module load hisat/2.1.0
  module load stringtie/1.3.3b
  module load xz/5.2.2
  module load htslib/1.3.3
  module load gffread/2
  module load perl/5.26.1
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The following functions are called from the reads2vcf.sh script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
source ./reads2vcf.sh

  # prepare environment
  loadModules
  createWorkingEnvironment-reads2vcf
  
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
  copyRef-reads2vcf
  
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
  
  ## -- rename files; see reads2vcf.sh for details
  renameSortedBAMs
 
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The following functions are called from the SNP_analysis.sh script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
source ./SNP_analysis.sh

  ## -- create full dataset (combining Jessica's WGS data with Seq-Cap+RNA-Seq)
    ### --- this function will create three files: Full_CDS.vcf, Full_Exons.vcf, 
    ### --- and Full_Genes.vcf
  combineDatasets
  
  # -- sort variants and remove duplicates
  sortVariants Full_CDS Full_IndividualsToRemove
  sortVariants Full_Exons Full_IndividualsToRemove
  sortVariants Full_Genes Full_IndividualsToRemove
  
  # -- annotate variants using snpeff
  functionalAnnotation Full
  
  # -- extract SNPs by gene
  getGeneVariants CDS
  getGeneVariants Exons
  getGeneVariants Genes
  getGeneVariants CDS _missense
  getGeneVariants CDS _synonymous
  
  # -- get transcript lengths and number of SNPs for each gene
  getTranscriptLengths CDS
  getTranscriptLengths Exons
  getTranscriptLengths Genes
  getTranscriptLengths CDS _missense
  getTranscriptLengths CDS _synonymous
  
  # -- convert vcf to fasta for targeted genes
  vcf2faa
  reference2faa
  # -- move captured target genes to new directory
  moveCapturedGenes
  # -- create alignments for peptide and nucleotide sequences
  createMSA faa protein Sequences maskedMSA
  createMSA fna transcript Sequences maskedMSA
  
  # -- create population text files for pops and pairwise comparisons
  createPopFiles
  # -- create vcf files containing samples for each pairwise pop comparison
  createPairwiseVCFs
  # -- calculate Tajima's D for each pairwise pop comparison for each gene
  getPairwisePopGen
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The following functions are called from the annotateVCF.sh script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
source ./sift2vcf.sh

  WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory/SNP_analysis/proteinStructure/Thamnophis_elegans
  # -- Test SIFT4G (optional- uncomment if you want to do it)
  # checkSIFT4G
  
  # -- create working environment
  createWorkingEnvironment-sift2vcf

  # -- copy reference genome and annotation
  copyRef-sift2vcf

  # -- download protein database
  downloadUniRef

  # ***** MAKE SURE THE CONFIGURATION FILE IS READY BEFORE CONTINUING TO NEXT STEP ***** #
  ## **********         see $WorkingDirectory/Thamnophis_elegans.txt        *********** ##

  # -- create genomic SIFT database using SIFT4G
  createSIFTdatabase
  
  # -- annotate vcf using SIFT database
  annotateVCF-sift2vcf
  # -- Annotate vcf with SIFT scores
