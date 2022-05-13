#!/bin/sh
#SBATCH --job-name=Thamnophis_main_tmp
#SBATCH --nodes=1                     # node(s) required for job
#SBATCH --ntasks=4                    # number of tasks across all nodes
#SBATCH --partition=jro0014_amd       # name of partition
#SBATCH --time=20:00:00               # Run time (D-HH:MM:SS)
#SBATCH --output=test-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=test_error-%j.err     # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL               # will send email for begin,end,fail
#SBATCH --mail-user=rlk0015@auburn.edu

WorkingDirectory=/scratch/rlk0015/Telegans_SeqCapProject/WorkingDirectory
pythonScripts=/home/rlk0015/projects/Thamnophis/ThamnophisElegans_FunctionalGenomics2021/Examining-Sequence-Variation/pythonScripts

function loadModules_Hopper {
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

function loadModules {
  module load bcftools/1.11
  module load samtools/1.11
  module load vcftools/0.1.17
  module load bedtools/2.30.0
  module load gffreader/12.7
  module load gatk/4.1.9.0
  conda init bash
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The following functions are called from the reads2vcf.sh script
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cd /home/rlk0015/projects/Thamnophis/ThamnophisElegans_FunctionalGenomics2021/Examining-Sequence-Variation
# source reads2vcf.sh

 # prepare environment
 # for Hopper: loadModules_Hopper
 loadModules # If you encounter errors, see if modules from loadModules_Hopper need inclusion
 # createWorkingEnvironment-reads2vcf
 
#+ COMPLETE  ## ------------------------
#+ COMPLETE  # begin raw read processing
#+ COMPLETE  ## ------------------------
#+ COMPLETE  # Copy Sequence Capture raw reads into working environment
#+ COMPLETE  copyRawReadsDNA
#+ COMPLETE  # Quality check RNA reads
#+ COMPLETE  performFASTQC rawReadsDNA
#+ COMPLETE  ## Note: RNA Seq reads were already copied
#+ COMPLETE  ## and quality checked. If this needs to be
#+ COMPLETE  ## done again, refer to the copyRawReadsDNA
#+ COMPLETE  ## and performFastQC commands in reads2vcf.sh
#+ COMPLETE  
#+ COMPLETE  # quality clean reads
#+ COMPLETE  performTrimmingPE DNA
#+ COMPLETE  performTrimmingPE RNA
#+ COMPLETE  performTrimmingSE RNA
#+ COMPLETE  
#+ COMPLETE  # Quality check clean reads
#+ COMPLETE  performFASTQC cleanReadsDNA
#+ COMPLETE  performFASTQC cleanReadsRNA
#+ COMPLETE  
#+ COMPLETE  # Copy reference genome into working environment
#+ COMPLETE  copyRef-reads2vcf
#+ COMPLETE  
#+ COMPLETE  # Map clean reads to reference genome
#+ COMPLETE  mapReadsDNA
#+ COMPLETE  mapPEReadsRNA
#+ COMPLETE  mapSEReadsRNA
#+ COMPLETE  
#+ COMPLETE  # Change the sequence capture names to match sample ID
#+ COMPLETE  changeSeqCapNames
#+ COMPLETE  
#+ COMPLETE  # Add sequencing read group information to each sample
#+ COMPLETE  AddReadGroupsRNA 
#+ COMPLETE  AddReadGroupsDNA
#+ COMPLETE  
#+ COMPLETE  # Prep for SNP calling and calculate mapping stats
#+ COMPLETE  indexReference TelegansGenome.fasta
#+ COMPLETE  prepForVariantCalling DNA
#+ COMPLETE  prepForVariantCalling RNA
#+ COMPLETE  
#+ COMPLETE  ### --------------------- Variant Calling ----------------------- ###
#+ COMPLETE  # Perform bqsr-'bootstraping', SNP calling, and hard filtration on Seq Cap data
#+ COMPLETE  cd $WorkingDirectory/GATKDNA
#+ COMPLETE  ## Replicate 1
#+ COMPLETE    use-HaplotypeCaller 0 DNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 0 
#+ COMPLETE    initial-VariantFiltration JustSNPs_0.vcf 0
#+ COMPLETE    use-BaseRecalibrator 0 DNA
#+ COMPLETE    use-BQSR 0 1 DNA
#+ COMPLETE    use-HaplotypeCaller 1 DNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 1
#+ COMPLETE    use-BaseRecalibrator 1 DNA
#+ COMPLETE    use-AnalyzeCovariates 0 1 DNA
#+ COMPLETE  ## -- Replicate 2
#+ COMPLETE    use-BQSR 1 2 DNA
#+ COMPLETE    use-HaplotypeCaller 2 DNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 2
#+ COMPLETE    use-BaseRecalibrator 2 DNA
#+ COMPLETE    use-AnalyzeCovariates 1 2 DNA
#+ COMPLETE  
#+ COMPLETE  # Perform bqsr-'bootstraping', SNP calling, and hard filtration on RNA-seq data
#+ COMPLETE  cd $WorkingDirectory/GATKRNA
#+ COMPLETE  ## Replicate 1
#+ COMPLETE    use-HaplotypeCaller 0 RNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 0
#+ COMPLETE    initial-VariantFiltration JustSNPs_0.vcf 0
#+ COMPLETE    use-BaseRecalibrator 0 RNA
#+ COMPLETE    use-BQSR 0 1 RNA
#+ COMPLETE    use-HaplotypeCaller 1 RNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 1
#+ COMPLETE    use-BaseRecalibrator 1 RNA
#+ COMPLETE    use-AnalyzeCovariates 0 1 RNA
#+ COMPLETE  ## -- Replicate 2
#+ COMPLETE    use-BQSR 1 2 RNA
#+ COMPLETE    use-HaplotypeCaller 2 RNA TelagGenome.fasta
#+ COMPLETE    get-just-SNPs 2
#+ COMPLETE    use-BaseRecalibrator 2 RNA
#+ COMPLETE    use-AnalyzeCovariates 1 2 RNA
#+ COMPLETE  
#+ COMPLETE  ## -- Merge RNA and DNA data
#+ COMPLETE  combine-VCF
#+ COMPLETE  cd $WorkingDirectory/variantFiltration 
#+ COMPLETE   ## -- Annotate variants
#+ COMPLETE   getNetworkFasta IILS
#+ COMPLETE   probes2gff Exons_2021.fa SeqCap
#+ COMPLETE   annotateVariants Merged.vcf.gz SeqCap
#+ COMPLETE   ## -- Initial Filter Variants
#+ COMPLETE   initial-VariantFiltration SeqCap_Annotated.vcf SeqCap_InitialFiltered
#+ COMPLETE   ## -- Perform Hard Filtering
#+ COMPLETE   hard-VariantFiltration SeqCap_InitialFiltered SeqCap
#+ COMPLETE   ## -- Eliminate potential RNA editing sites
#+ COMPLETE   removeRNAedits
#+ COMPLETE   ## -- Get coding regions and exonic regions
#+ COMPLETE   getSpecificVariants SeqCap CDS
#+ COMPLETE   getSpecificVariants SeqCap Exons
#+ COMPLETE   ## -- Copy final step of hard filtering (for nomenclature purposes)
#+ COMPLETE   cp SeqCap_HardFilterStep6.vcf SeqCap_Genes.vcf
#+ COMPLETE   ## -- plot the variants (you can do this with all of the files as desired)
#+ COMPLETE   plotVariants SeqCap_Exons.vcf
#+ COMPLETE   plotVariants SeqCap_Annotated.vcf
#+ COMPLETE   ## -- zip up everything
#+ COMPLETE   bgzip SeqCap_Genes.vcf
#+ COMPLETE   bgzip SeqCap_Exons.vcf
#+ COMPLETE   bgzip SeqCap_CDS.vcf
#+ COMPLETE    
#+ COMPLETE  ## -- rename files; see reads2vcf.sh for details
#+ COMPLETE  renameSortedBAMs
#+ COMPLETE  
#+ COMPLETE   
#+ COMPLETE  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#+ COMPLETE  # The following functions are called from the SNP_curation.sh script
#+ COMPLETE  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cd /home/rlk0015/projects/Thamnophis/ThamnophisElegans_FunctionalGenomics2021/Examining-Sequence-Variation
source ./SNP_curation.sh
#+ COMPLETE  
#+ COMPLETE    cd $WorkingDirectory/variantFiltration
#+ COMPLETE    ## -- create full dataset (combining Jessica's WGS data with Seq-Cap+RNA-Seq)
#+ COMPLETE      ### --- this function will create three files: Full_CDS.vcf, Full_Exons.vcf, 
#+ COMPLETE      ### --- and Full_Genes.vcf
#+ COMPLETE    combineDatasets
#+ COMPLETE    
#+ COMPLETE    # -- sort variants and remove duplicates
#+ COMPLETE    sortVariants Full_CDS $WorkingDirectory/References/Full_IndividualsToRemove.txt
#+ COMPLETE    sortVariants Full_Exons $WorkingDirectory/References/Full_IndividualsToRemove.txt
#+ COMPLETE    sortVariants Full_Genes $WorkingDirectory/References/Full_IndividualsToRemove.txt
#+ COMPLETE    
#+ COMPLETE  # -- annotate variants using snpeff
#+ COMPLETE  functionalAnnotation Full
#+ COMPLETE  
#+ COMPLETE  # -- extract SNPs by gene
#+ COMPLETE  getGeneVariants CDS
#+ COMPLETE  getGeneVariants Exons
#+ COMPLETE  getGeneVariants Genes
#+ COMPLETE  getGeneVariants CDS _missense
#+ COMPLETE  getGeneVariants CDS _synonymous
#+ COMPLETE  
#+ COMPLETE  # -- get transcript lengths and number of SNPs for each gene
#+ COMPLETE  getTranscriptLengths CDS
#+ COMPLETE  getTranscriptLengths Exons
#+ COMPLETE  getTranscriptLengths Genes
#+ COMPLETE  getTranscriptLengths CDS _missense
#+ COMPLETE  getTranscriptLengths CDS _synonymous
#+ COMPLETE  
#+ COMPLETE  # -- convert vcf to fasta for targeted genes
#+ COMPLETE  vcf2faa
#+ COMPLETE  reference2faa
#+ COMPLETE  # -- move captured target genes to new directory
#+ COMPLETE  moveCapturedGenes
#+ COMPLETE  # -- create alignments for peptide and nucleotide sequences
#+ COMPLETE  createMSA faa protein Sequences maskedMSA
#+ COMPLETE  createMSA fna transcript Sequences maskedMSA
#+ COMPLETE 
#+ COMPLETE  # -- create population text files for pops and pairwise comparisons
#+ COMPLETE  createPopFiles
#+ COMPLETE  # -- create vcf files containing samples for each pairwise pop comparison
#+ COMPLETE  createPairwiseVCFs
#+ COMPLETE  getGeneTajD Exons
#+ COMPLETE  pairwisePopGen2
#+ COMPLETE  sort_samples_by_ecotype Full_Genes
#+ COMPLETE  get_ecotype_vcf Full_Genes Lakeshore 
#+ COMPLETE  get_ecotype_vcf Full_Genes Meadow 
#+ COMPLETE  get_ecotype_allele_freq Full_Genes
#+ COMPLETE  sort_samples_by_ecotype top_fst_dxy_sift
#+ COMPLETE  get_ecotype_vcf top_fst_dxy_sift Lakeshore 
#+ COMPLETE  get_ecotype_vcf top_fst_dxy_sift Meadow 
#+ COMPLETE  get_ecotype_allele_freq top_fst_dxy_sift
# -- calculate PopGen statistics for each pairwise pop comparison
#+ NO LONGER DOING: getPairwisePopGen

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
