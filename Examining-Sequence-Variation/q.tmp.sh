#!/bin/sh

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
cd /home/rlk0015/SeqCap/code/GenomicProcessingPipeline/Examining-Sequence-Variation
source reads2vcf.sh

  # prepare environment
  loadModules
#+ COMPLETED   createWorkingEnvironment-reads2vcf
#+ COMPLETED   
#+ COMPLETED   ## ------------------------
#+ COMPLETED   # begin raw read processing
#+ COMPLETED   ## ------------------------
#+ COMPLETED   # Copy Sequence Capture raw reads into working environment
#+ COMPLETED   copyRawReadsDNA
#+ COMPLETED   # Quality check RNA reads
#+ COMPLETED   performFASTQC rawReadsDNA
#+ COMPLETED   ## Note: RNA Seq reads were already copied
#+ COMPLETED   ## and quality checked. If this needs to be
#+ COMPLETED   ## done again, refer to the copyRawReadsDNA
#+ COMPLETED   ## and performFastQC commands in reads2vcf.sh
#+ COMPLETED   
#+ COMPLETED   # quality clean reads
#+ COMPLETED   performTrimmingPE DNA
#+ COMPLETED   performTrimmingPE RNA
#+ COMPLETED   performTrimmingSE RNA
#+ COMPLETED   
#+ COMPLETED   # Quality check clean reads
#+ COMPLETED   performFASTQC cleanReadsDNA
#+ COMPLETED   performFASTQC cleanReadsRNA
#+ COMPLETED   
#+ COMPLETED   # Copy reference genome into working environment
#+ COMPLETED   copyRef-reads2vcf
#+ COMPLETED   
#+ COMPLETED   # Map clean reads to reference genome
#+ COMPLETED   mapReadsDNA
#+ COMPLETED   mapPEReadsRNA
#+ COMPLETED   mapSEReadsRNA
#+ COMPLETED   
#+ COMPLETED   # Change the sequence capture names to match sample ID
#+ COMPLETED   changeSeqCapNames
#+ COMPLETED   
#+ COMPLETED   # Add sequencing read group information to each sample
#+ COMPLETED   readGroupsRNA 
#+ COMPLETED   readGroupsDNA
#+ COMPLETED   
#+ COMPLETED   # Prep for SNP calling and calculate mapping stats
#+ COMPLETED   indexReference
#+ COMPLETED   prepForVariantCalling DNA
#+ COMPLETED   prepForVariantCalling RNA
#+ COMPLETED   
#+ COMPLETED   ### --------------------- Variant Calling ----------------------- ###
#+ COMPLETED   # Perform bqsr-'bootstraping', SNP calling, and hard filtration on Seq Cap data
#+ COMPLETED   cd $WorkingDirectory/GATKDNA
#+ COMPLETED   ## Replicate 1
#+ COMPLETED     use-HaplotypeCaller 0 DNA
#+ COMPLETED     get-just-SNPs 0 
#+ COMPLETED     initial-VariantFiltration JustSNPs_0.vcf 0
#+ COMPLETED     use-BaseRecalibrator 0 DNA
#+ COMPLETED     use-BQSR 0 1 DNA
#+ COMPLETED     use-HaplotypeCaller 1 DNA
#+ COMPLETED     get-just-SNPs 1
#+ COMPLETED     use-BaseRecalibrator 1 DNA
#+ COMPLETED     use-AnalyzeCovariates 0 1 DNA
#+ COMPLETED   ## -- Replicate 2
#+ COMPLETED     use-BQSR 1 2 DNA
#+ COMPLETED     use-HaplotypeCaller 2 DNA
#+ COMPLETED     get-just-SNPs 2
#+ COMPLETED     use-BaseRecalibrator 2 DNA
#+ COMPLETED     use-AnalyzeCovariates 1 2 DNA
#+ COMPLETED   
#+ COMPLETED   # Perform bqsr-'bootstraping', SNP calling, and hard filtration on RNA-seq data
#+ COMPLETED   cd $WorkingDirectory/GATKRNA
#+ COMPLETED   ## Replicate 1
#+ COMPLETED     use-HaplotypeCaller 0 RNA
#+ COMPLETED     get-just-SNPs 0
#+ COMPLETED     initial-VariantFiltration JustSNPs_0.vcf 0
#+ COMPLETED     use-BaseRecalibrator 0 RNA
#+ COMPLETED     use-BQSR 0 1 RNA
#+ COMPLETED     use-HaplotypeCaller 1 RNA
#+ COMPLETED     get-just-SNPs 1
#+ COMPLETED     use-BaseRecalibrator 1 RNA
#+ COMPLETED     use-AnalyzeCovariates 0 1 RNA
#+ COMPLETED   ## -- Replicate 2
#+ COMPLETED     use-BQSR 1 2 RNA
#+ COMPLETED     use-HaplotypeCaller 2 RNA
#+ COMPLETED     get-just-SNPs 2
#+ COMPLETED     use-BaseRecalibrator 2 RNA
#+ COMPLETED     use-AnalyzeCovariates 1 2 RNA
#+ COMPLETED   
#+ COMPLETED   ## -- Merge RNA and DNA data
#+ COMPLETED   combine-VCF
#+ COMPLETED   cd $WorkingDirectory/variantFiltration 
#+ COMPLETED   ## -- Eliminate potential RNA editing sites
  cd $WorkingDirectory/variantFiltration
  removeRNAedits
  ## -- Annotate variants
  #+ COMPLETED getNetworkFasta IILS
  #+ COMPLETED probes2gff Exons_2021.fa SeqCap
  annotateVariants removedRNAedits.vcf SeqCap
  ## -- Initial Filter Variants
  initial-VariantFiltration SeqCap_Annotated.vcf SeqCap_InitialFiltered
  ## -- Perform Hard Filtering
  hard-VariantFiltration SeqCap_InitialFiltered SeqCap
  ## -- Get coding regions and exonic regions
  getSpecificVariants SeqCap CDS
  getSpecificVariants SeqCap Exons
  ## -- Copy final step of hard filtering (for nomenclature purposes)
  cp SeqCap_HardFilterStep4.vcf SeqCap_Genes.vcf
  ## -- zip up everything
  bgzip SeqCap_Genes.vcf
  bgzip SeqCap_Exons.vcf
  bgzip SeqCap_CDS.vcf
  ## -- plot the variants (you can do this with all of the files as desired)
  plotVariants SeqCap_Exons.vcf
  plotVariants SeqCap_Annotated.vcf
  
  #+ COMPLETED ## -- rename files; see reads2vcf.sh for details
  #+ COMPLETED renameSortedBAMs
 
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The following functions are called from the SNP_curation.sh script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cd /home/rlk0015/SeqCap/code/GenomicProcessingPipeline/Examining-Sequence-Variation
source SNP_curation.sh

  cd $WorkingDirectory/variantFiltration
  ## -- create full dataset (combining Jessica's WGS data with Seq-Cap+RNA-Seq)
    ### --- this function will create three files: Full_CDS.vcf, Full_Exons.vcf, 
    ### --- and Full_Genes.vcf
  combineDatasets
  
  # -- sort variants and remove duplicates
  sortVariants Full_CDS $WorkingDirectory/Full_IndividualsToRemove.txt
  sortVariants Full_Exons $WorkingDirectory/Full_IndividualsToRemove.txt
  sortVariants Full_Genes $WorkingDirectory/Full_IndividualsToRemove.txt
  
#+ WAITING  # -- annotate variants using snpeff
#+ WAITING  functionalAnnotation Full
#+ WAITING  
#+ WAITING  # -- extract SNPs by gene
#+ WAITING  getGeneVariants CDS
#+ WAITING  getGeneVariants Exons
#+ WAITING  getGeneVariants Genes
#+ WAITING  getGeneVariants CDS _missense
#+ WAITING  getGeneVariants CDS _synonymous
#+ WAITING  
#+ WAITING  # -- get transcript lengths and number of SNPs for each gene
#+ WAITING  getTranscriptLengths CDS
#+ WAITING  getTranscriptLengths Exons
#+ WAITING  getTranscriptLengths Genes
#+ WAITING  getTranscriptLengths CDS _missense
#+ WAITING  getTranscriptLengths CDS _synonymous
#+ WAITING  
#+ WAITING  # -- convert vcf to fasta for targeted genes
#+ WAITING  vcf2faa
#+ WAITING  reference2faa
#+ WAITING  # -- move captured target genes to new directory
#+ WAITING  moveCapturedGenes
#+ WAITING  # -- create alignments for peptide and nucleotide sequences
#+ WAITING  createMSA faa protein Sequences maskedMSA
#+ WAITING  createMSA fna transcript Sequences maskedMSA
#+ WAITING  
#+ WAITING  # -- create population text files for pops and pairwise comparisons
#+ WAITING  createPopFiles
#+ WAITING  # -- create vcf files containing samples for each pairwise pop comparison
#+ WAITING  createPairwiseVCFs
#+ WAITING  # -- calculate Tajima's D for each pairwise pop comparison for each gene
#+ WAITING  getPairwisePopGen
#+ WAITING  
#+ WAITING# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#+ WAITING# The following functions are called from the annotateVCF.sh script
#+ WAITING# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#+ WAITINGsource ./sift2vcf.sh
#+ WAITING
#+ WAITING  WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory/SNP_analysis/proteinStructure/Thamnophis_elegans
#+ WAITING  # -- Test SIFT4G (optional- uncomment if you want to do it)
#+ WAITING  # checkSIFT4G
#+ WAITING  
#+ WAITING  # -- create working environment
#+ WAITING  createWorkingEnvironment-sift2vcf
#+ WAITING
#+ WAITING  # -- copy reference genome and annotation
#+ WAITING  copyRef-sift2vcf
#+ WAITING
#+ WAITING  # -- download protein database
#+ WAITING  downloadUniRef
#+ WAITING
#+ WAITING  # ***** MAKE SURE THE CONFIGURATION FILE IS READY BEFORE CONTINUING TO NEXT STEP ***** #
#+ WAITING  ## **********         see $WorkingDirectory/Thamnophis_elegans.txt        *********** ##
#+ WAITING
#+ WAITING  # -- create genomic SIFT database using SIFT4G
#+ WAITING  createSIFTdatabase
#+ WAITING  
#+ WAITING  # -- annotate vcf using SIFT database
#+ WAITING  annotateVCF-sift2vcf
#+ WAITING  # -- Annotate vcf with SIFT scores
