#!/bin/sh

#Give job a name
#PBS -N FullScript_Feb2019

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=2:ppn=8,walltime=40:00:00 
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

module load fastqc/11.5
module load gnu-parallel/20160322 
module load trimmomatic/0.36
module load bwa/0.7.15
module load samtools/1.3.1
module load picard/2.4.1
module load python/2.7.15
module load gatk/3.6
module load bcftools/1.3.2
module load hybpiper/1
module load xz/5.2.2
module load htslib/1.3.3

#create working directory ###
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory
# DNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/assembledReads
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA
# RNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsRNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
echo "RLK_report: directory created: /scratch/rlk0015/Telag/Dec2018/WorkingDirectory with rawReads and cleanReads sub directories"

#### The following is commented out because I already have cleaned reads in the HybPiperContigs working directory. Remove comment when publishing ####
# -----------------------------------------------------------------------------------
# ### copy raw data to scratch ###
# cd /scratch/rlk0015/Telag/rawReadsDNA
# ls HTAdapter*.fastq.gz | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA'
# 
# ### quality check ###
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA
# ls *.fastq | time parallel -j+0 --eta 'fastqc {}'
# 
# ### copy over adapter file ###
# cp /home/rlk0015/SeqCap/code/References/adapters.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA
# 
# ### paired-end trimming ###
# ls | grep "fastq.gz" | cut -d "_" -f 1,2 | sort | uniq > PE_TrimmList
# while read i
# do
# java -jar /tools/trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA/"$i"*_R1*.fastq.gz /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/rawReadsDNA/"$i"*_R2*.fastq.gz /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_R1_paired.fastq.gz /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_R1_unpaired.fastq.gz /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_R2_paired.fastq.gz /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_R2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:30 MINLEN:36
# done<PE_TrimmList
# ### quality check ###
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/
# ls *paired.fastq.gz | grep "paired.fastq.gz" | time parallel -j+0 --eta 'fastqc {}'
# -----------------------------------------------------------------------------------

#### ************** COPY STUFF OVER ************** ####
# Copy cleaned reads over
# ---------------------------
# DONE:
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory_HybPiperContigs_T2/cleanReadsDNA
# ls *.fastq | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA'
# ---------------------------
# Copy clean RNA reads over
# ---------------------------
# Don't need for now (already mapped)
# cd /scratch/GarterSnake/RNAseq_2012Data/CleanedData
# ls SRR*.fastq | parallel -j+0 --eta 'cp {} /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA'
# mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/Report/RNAseq_2012copied
# cd /scratch/GarterSnake/RNAseq_2008Data/CleanedData
# ls SRR*.fastq | parallel -j+0 --eta 'cp {} /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA'
# mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/Report/RNAseq_2008copied
# # Copy mapped reads over
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory_HybPiperContigs_T2/mappedReadsRNA
# ls *.sam | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA'
# mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/Report/mappedReadsRNA_copied
#### ********************************************* ####

### -------- HybPiper Assembly -------- ###
# ----------------------------------
# This is already done:
# cp /home/rlk0015/SeqCap/code/References/Transcripts.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/assembledReads
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA
# ls | grep "fastq" | cut -d "_" -f 1 | sort | uniq > cleanReadsList
# while read i
# do
# reads_first.py -b /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/assembledReads/Transcripts.fa -r /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"*_R*_paired.fastq --bwa --cpu 20 --prefix /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/assembledReads/HybPiperAssembly
# done<cleanReadsList
# # Create test_seq_lengths.txt file (for data visualization)
# echo "HybPiperAssembly" >> names.txt
# python /tools/hybpiper-1/get_seq_lengths.py Transcripts.fa names.txt dna > test_seq_lengths.txt
# # Retrieve Assembled Sequences
# python /tools/hybpiper-1/retrieve_sequences.py Transcripts.fa . dna
# mkdir -p Results
# mv *.FNA Results
# cd Results
# ls *.FNA | cut -d "." -f 1 | sort > list
# while read i
# do
# sed -i "s/>.\+/>Telegans-$i/" "$i".FNA
# done<list
# # Create reference
# cat *.FNA > HybPiper_Transcripts_Round1.fasta
# cp HybPiper_Transcripts_Round1.fasta /home/rlk0015/SeqCap/code/References
# ----------------------------------

# ------------------------------------------------------------------
# Tonia's criteria for trimming RNA-Seq reads:
############ 2008 Data Trimmomatic #############
# java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar SE -phred33 "$i"_1.fastq  /scratch/GarterSnake/RNAseq_2008Data/CleanedData/"$i"_cleaned.fastq  LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:36 2012 samples are PE 100 bp reads.
############ 2012 Trimmomatic #############
# java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar PE  -phred33 "$i"_1.fastq "$i"_2.fastq "$i"_1_paired.fastq "$i"_1_unpaired.fastq "$i"_2_paired.fastq "$i"_2_unpaired.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
# ------------------------------------------------------------------

# ******************************* Completed the following, uncommment for publication ************************************
#  # Quality check on RNA-Seq data
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA
#  ls *.fastq | grep "fastq" | time parallel -j+0 --eta 'fastqc {}'
#  
#  ### ***  -------------------  MAP DNA TO REFERENCE  ------------------ *** ###
#  mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper
#  mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts
#  ### copy the reference to mapped reads directory ###
#  cp /home/rlk0015/SeqCap/code/References/HybPiper_Transcripts_Round1.fasta /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper/HybPiperContigs.fa
#  cp /home/rlk0015/SeqCap/code/References/Transcripts.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts/Transcripts.fa
#  echo "RLK_report: reference copy complete"
#  # index reference
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper
#  bwa index -p hybPiperContigs -a is HybPiperContigs.fa
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts
#  bwa index -p transcripts -a is Transcripts.fa
#  echo "RLK_report: reference index complete"
#  # create list with each paired individual
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA
#  ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/pairedMapList
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA
#  # while loop through the names in pairedMapList
#  while read i
#  do
#  ### map to ref hybPiperContigs ###
#  bwa mem -t 4 -M HybPiper/hybPiperContigs /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_*R1_paired.fastq /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_*R2_paired.fastq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper/"$i"_hybPiper_doubles.sam
#  bwa mem -t 4 -M Transcripts/transcripts /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_*R1_paired.fastq /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsDNA/"$i"_*R2_paired.fastq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"_transcripts_doubles.sam
#  done<pairedMapList
#  echo "RLK_report: map complete"
#  
#  # ********************* THIS IS COMPLETE (uncomment when publishing) ***************************
#  # ### ***  -------------------  MAP RNA TO REFERENCE  ------------------ *** ###
#  # # ********* Already done ***********
#  # ### copy the reference to mapped reads directory ###
#  # cp /home/rlk0015/SeqCap/code/References/Transcriptome.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/Transcriptome.fa
#  # echo "RLK_report: reference copy complete"
#  # # index reference
#  # cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
#  # bwa index -p transcriptome -a is Transcriptome.fa
#  # echo "RLK_report: reference index complete"
#  # # **********************************
#  # 
#  # # create list with each paired individual
#  # cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA
#  # ls | grep "_cleaned.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/singleEndMapList
#  # cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
#  # # ----------------------------------
#  # # Copied over mapped reads, so below is unnecessary for now
#  # # # while loop through the names in pairedEndMapList
#  # # while read i
#  # # do
#  # # ### PAIRED-END MAPPING ###
#  # # bwa mem -t 4 -M transcriptome /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA/"$i"*_1_paired.fastq $# WorkingDirectory/cleanReadsRNA/"$i"*_2_paired.fastq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory# /mappedReadsRNA/"$i"_doubles.sam
#  # # done</scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/pairedEndMapList
#  # # echo "RLK_report: paired map complete"
#  # # ----------------------------------
#  # ### SINGLE-END MAPPING ###
#  # while read i
#  # do
#  # bwa mem -t 4 -M transcriptome /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/cleanReadsRNA/"$i"*.fastq > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/"$i"_doubles.sam
#  # done</scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/singleEndMapList
#  # 
#  # ***********************************************************************
#  
#  ### change names to match population sampling ###
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper
#  mv HTAdapter1_hybPiper_doubles.sam ELF01_hybPiper_doubles.sam
#  mv HTAdapter2_hybPiper_doubles.sam ELF02_hybPiper_doubles.sam
#  mv HTAdapter3_hybPiper_doubles.sam ELF03_hybPiper_doubles.sam
#  mv HTAdapter4_hybPiper_doubles.sam ELF04_hybPiper_doubles.sam
#  mv HTAdapter5_hybPiper_doubles.sam ELF05_hybPiper_doubles.sam
#  mv HTAdapter6_hybPiper_doubles.sam ELF06_hybPiper_doubles.sam
#  mv HTAdapter7_hybPiper_doubles.sam ELF07_hybPiper_doubles.sam
#  mv HTAdapter8_hybPiper_doubles.sam ELF08_hybPiper_doubles.sam
#  mv HTAdapter9_hybPiper_doubles.sam ELF09_hybPiper_doubles.sam
#  mv HTAdapter10_hybPiper_doubles.sam ELF10_hybPiper_doubles.sam
#  mv HTAdapter11_hybPiper_doubles.sam ELF11_hybPiper_doubles.sam
#  mv HTAdapter12_hybPiper_doubles.sam ELF12_hybPiper_doubles.sam
#  mv HTAdapter13_hybPiper_doubles.sam ELF13_hybPiper_doubles.sam
#  mv HTAdapter14_hybPiper_doubles.sam ELF14_hybPiper_doubles.sam
#  mv HTAdapter15_hybPiper_doubles.sam ELF15_hybPiper_doubles.sam
#  mv HTAdapter16_hybPiper_doubles.sam ELF16_hybPiper_doubles.sam
#  mv HTAdapter17_hybPiper_doubles.sam MAH01_hybPiper_doubles.sam
#  mv HTAdapter18_hybPiper_doubles.sam MAH02_hybPiper_doubles.sam
#  mv HTAdapter19_hybPiper_doubles.sam MAH03_hybPiper_doubles.sam
#  mv HTAdapter20_hybPiper_doubles.sam MAH04_hybPiper_doubles.sam
#  mv HTAdapter21_hybPiper_doubles.sam MAH05_hybPiper_doubles.sam
#  mv HTAdapter22_hybPiper_doubles.sam MAH06_hybPiper_doubles.sam
#  mv HTAdapter23_hybPiper_doubles.sam MAH07_hybPiper_doubles.sam
#  mv HTAdapter24_hybPiper_doubles.sam MAH08_hybPiper_doubles.sam
#  mv HTAdapter25_hybPiper_doubles.sam MAH09_hybPiper_doubles.sam
#  mv HTAdapter26_hybPiper_doubles.sam MAH10_hybPiper_doubles.sam
#  mv HTAdapter27_hybPiper_doubles.sam MAH11_hybPiper_doubles.sam
#  mv HTAdapter28_hybPiper_doubles.sam MAH12_hybPiper_doubles.sam
#  mv HTAdapter29_hybPiper_doubles.sam MAH13_hybPiper_doubles.sam
#  mv HTAdapter30_hybPiper_doubles.sam MAH14_hybPiper_doubles.sam
#  mv HTAdapter31_hybPiper_doubles.sam MAH15_hybPiper_doubles.sam
#  mv HTAdapter32_hybPiper_doubles.sam MAR01_hybPiper_doubles.sam
#  mv HTAdapter33_hybPiper_doubles.sam MAR02_hybPiper_doubles.sam
#  mv HTAdapter34_hybPiper_doubles.sam MAR03_hybPiper_doubles.sam
#  mv HTAdapter35_hybPiper_doubles.sam MAR04_hybPiper_doubles.sam
#  mv HTAdapter36_hybPiper_doubles.sam MAR05_hybPiper_doubles.sam
#  mv HTAdapter37_hybPiper_doubles.sam MAR06_hybPiper_doubles.sam
#  mv HTAdapter38_hybPiper_doubles.sam MAR07_hybPiper_doubles.sam
#  mv HTAdapter39_hybPiper_doubles.sam MAR08_hybPiper_doubles.sam
#  mv HTAdapter40_hybPiper_doubles.sam MAR09_hybPiper_doubles.sam
#  mv HTAdapter41_hybPiper_doubles.sam MAR10_hybPiper_doubles.sam
#  mv HTAdapter42_hybPiper_doubles.sam PAP01_hybPiper_doubles.sam
#  mv HTAdapter43_hybPiper_doubles.sam PAP02_hybPiper_doubles.sam
#  mv HTAdapter44_hybPiper_doubles.sam PAP03_hybPiper_doubles.sam
#  mv HTAdapter45_hybPiper_doubles.sam PAP04_hybPiper_doubles.sam
#  mv HTAdapter46_hybPiper_doubles.sam PAP05_hybPiper_doubles.sam
#  mv HTAdapter47_hybPiper_doubles.sam PAP06_hybPiper_doubles.sam
#  mv HTAdapter48_hybPiper_doubles.sam PAP07_hybPiper_doubles.sam
#  mv HTAdapter49_hybPiper_doubles.sam PAP08_hybPiper_doubles.sam
#  mv HTAdapter50_hybPiper_doubles.sam PAP09_hybPiper_doubles.sam
#  mv HTAdapter51_hybPiper_doubles.sam PAP10_hybPiper_doubles.sam
#  mv HTAdapter52_hybPiper_doubles.sam PAP11_hybPiper_doubles.sam
#  mv HTAdapter53_hybPiper_doubles.sam PAP12_hybPiper_doubles.sam
#  mv HTAdapter54_hybPiper_doubles.sam PAP13_hybPiper_doubles.sam
#  mv HTAdapter55_hybPiper_doubles.sam PAP14_hybPiper_doubles.sam
#  mv HTAdapter56_hybPiper_doubles.sam PAP15_hybPiper_doubles.sam
#  mv HTAdapter57_hybPiper_doubles.sam STO01_hybPiper_doubles.sam
#  mv HTAdapter58_hybPiper_doubles.sam STO02_hybPiper_doubles.sam
#  mv HTAdapter59_hybPiper_doubles.sam STO03_hybPiper_doubles.sam
#  mv HTAdapter60_hybPiper_doubles.sam STO04_hybPiper_doubles.sam
#  mv HTAdapter61_hybPiper_doubles.sam STO05_hybPiper_doubles.sam
#  mv HTAdapter62_hybPiper_doubles.sam STO06_hybPiper_doubles.sam
#  mv HTAdapter63_hybPiper_doubles.sam STO07_hybPiper_doubles.sam
#  mv HTAdapter64_hybPiper_doubles.sam STO08_hybPiper_doubles.sam
#  mv HTAdapter65_hybPiper_doubles.sam STO09_hybPiper_doubles.sam
#  mv HTAdapter66_hybPiper_doubles.sam STO10_hybPiper_doubles.sam
#  mv HTAdapter67_hybPiper_doubles.sam STO11_hybPiper_doubles.sam
#  mv HTAdapter68_hybPiper_doubles.sam STO12_hybPiper_doubles.sam
#  mv HTAdapter69_hybPiper_doubles.sam STO13_hybPiper_doubles.sam
#  mv HTAdapter70_hybPiper_doubles.sam STO14_hybPiper_doubles.sam
#  mv HTAdapter71_hybPiper_doubles.sam STO15_hybPiper_doubles.sam
#  mv HTAdapter72_hybPiper_doubles.sam STO17_hybPiper_doubles.sam
#  mv HTAdapter73_hybPiper_doubles.sam STO18_hybPiper_doubles.sam
#  mv HTAdapter74_hybPiper_doubles.sam STO19_hybPiper_doubles.sam
#  mv HTAdapter75_hybPiper_doubles.sam STO20_hybPiper_doubles.sam
#  mv HTAdapter76_hybPiper_doubles.sam STO21_hybPiper_doubles.sam
#  mv HTAdapter77_hybPiper_doubles.sam SUM01_hybPiper_doubles.sam
#  mv HTAdapter78_hybPiper_doubles.sam SUM02_hybPiper_doubles.sam
#  mv HTAdapter79_hybPiper_doubles.sam SUM03_hybPiper_doubles.sam
#  mv HTAdapter80_hybPiper_doubles.sam SUM04_hybPiper_doubles.sam
#  mv HTAdapter81_hybPiper_doubles.sam SUM05_hybPiper_doubles.sam
#  mv HTAdapter82_hybPiper_doubles.sam SUM06_hybPiper_doubles.sam
#  mv HTAdapter83_hybPiper_doubles.sam SUM07_hybPiper_doubles.sam
#  mv HTAdapter84_hybPiper_doubles.sam SUM08_hybPiper_doubles.sam
#  mv HTAdapter85_hybPiper_doubles.sam SUM09_hybPiper_doubles.sam
#  mv HTAdapter86_hybPiper_doubles.sam SUM10_hybPiper_doubles.sam
#  mv HTAdapter87_hybPiper_doubles.sam SUM11_hybPiper_doubles.sam
#  mv HTAdapter88_hybPiper_doubles.sam SUM12_hybPiper_doubles.sam
#  mv HTAdapter89_hybPiper_doubles.sam SUM13_hybPiper_doubles.sam
#  mv HTAdapter90_hybPiper_doubles.sam SUM15_hybPiper_doubles.sam
#  mv HTAdapter91_hybPiper_doubles.sam SUM16_hybPiper_doubles.sam
#  mv HTAdapter92_hybPiper_doubles.sam SUM17_hybPiper_doubles.sam
#  mv HTAdapter93_hybPiper_doubles.sam SUM18_hybPiper_doubles.sam
#  mv HTAdapter94_hybPiper_doubles.sam SUM19_hybPiper_doubles.sam
#  mv HTAdapter95_hybPiper_doubles.sam SUM20_hybPiper_doubles.sam
#  mv HTAdapter96_hybPiper_doubles.sam SUM21_hybPiper_doubles.sam
#  
#  ## -- ADD READGROUPS -- ##
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_hybPiper_doubles.sam" O="ELF01_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_54-38" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_hybPiper_doubles.sam" O="ELF02_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RA607" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_hybPiper_doubles.sam" O="ELF03_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_hybPiper_doubles.sam" O="ELF04_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_hybPiper_doubles.sam" O="ELF05_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_558-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_hybPiper_doubles.sam" O="ELF06_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_565-13" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_hybPiper_doubles.sam" O="ELF07_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_5651-15" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_hybPiper_doubles.sam" O="ELF08_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_568-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_hybPiper_doubles.sam" O="ELF09_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_546-03" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_hybPiper_doubles.sam" O="ELF10_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_547-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_hybPiper_doubles.sam" O="ELF11_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_6050-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_hybPiper_doubles.sam" O="ELF12_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_622-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_hybPiper_doubles.sam" O="ELF13_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_636-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_hybPiper_doubles.sam" O="ELF14_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_hybPiper_doubles.sam" O="ELF15_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_hybPiper_doubles.sam" O="ELF16_hybPiper_IDed.sam" RGPU="ELF" RGSM="ELF_RP54" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_hybPiper_doubles.sam" O="MAH01_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_35-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_hybPiper_doubles.sam" O="MAH02_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_520-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_hybPiper_doubles.sam" O="MAH03_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_539-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_hybPiper_doubles.sam" O="MAH04_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_540-07" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_hybPiper_doubles.sam" O="MAH05_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_541-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_hybPiper_doubles.sam" O="MAH06_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_542-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_hybPiper_doubles.sam" O="MAH07_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_543-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_hybPiper_doubles.sam" O="MAH08_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RA620" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_hybPiper_doubles.sam" O="MAH09_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RA641" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_hybPiper_doubles.sam" O="MAH10_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_hybPiper_doubles.sam" O="MAH11_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_hybPiper_doubles.sam" O="MAH12_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_hybPiper_doubles.sam" O="MAH13_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_hybPiper_doubles.sam" O="MAH14_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_hybPiper_doubles.sam" O="MAH15_hybPiper_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_hybPiper_doubles.sam" O="MAR01_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RP6771" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_hybPiper_doubles.sam" O="MAR02_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RP686" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_hybPiper_doubles.sam" O="MAR03_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_368- 2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_hybPiper_doubles.sam" O="MAR04_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_hybPiper_doubles.sam" O="MAR05_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA379" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_hybPiper_doubles.sam" O="MAR06_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_hybPiper_doubles.sam" O="MAR07_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_hybPiper_doubles.sam" O="MAR08_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA605" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_hybPiper_doubles.sam" O="MAR09_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA630" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_hybPiper_doubles.sam" O="MAR10_hybPiper_IDed.sam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_hybPiper_doubles.sam" O="PAP01_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_561-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_hybPiper_doubles.sam" O="PAP02_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_562-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_hybPiper_doubles.sam" O="PAP03_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_564-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_hybPiper_doubles.sam" O="PAP04_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_569-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_hybPiper_doubles.sam" O="PAP05_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_5010-05" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_hybPiper_doubles.sam" O="PAP06_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_5021-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_hybPiper_doubles.sam" O="PAP07_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_505-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_hybPiper_doubles.sam" O="PAP08_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_516-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_hybPiper_doubles.sam" O="PAP09_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_5241-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_hybPiper_doubles.sam" O="PAP10_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_28-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_hybPiper_doubles.sam" O="PAP11_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_RA434" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_hybPiper_doubles.sam" O="PAP12_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_RA647" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_hybPiper_doubles.sam" O="PAP13_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_RA648" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_hybPiper_doubles.sam" O="PAP14_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_RA649" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_hybPiper_doubles.sam" O="PAP15_hybPiper_IDed.sam" RGPU="PAP" RGSM="PAP_RP16" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_hybPiper_doubles.sam" O="STO01_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_RP697" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_hybPiper_doubles.sam" O="STO02_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_RP698" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_hybPiper_doubles.sam" O="STO03_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_RA530" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_hybPiper_doubles.sam" O="STO04_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_RA549" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_hybPiper_doubles.sam" O="STO05_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_RP22" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_hybPiper_doubles.sam" O="STO06_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_271" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_hybPiper_doubles.sam" O="STO07_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_275" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_hybPiper_doubles.sam" O="STO08_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-113" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_hybPiper_doubles.sam" O="STO09_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-115" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_hybPiper_doubles.sam" O="STO10_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-125" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_hybPiper_doubles.sam" O="STO11_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-76" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_hybPiper_doubles.sam" O="STO12_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-77" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_hybPiper_doubles.sam" O="STO13_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-78" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_hybPiper_doubles.sam" O="STO14_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-74" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_hybPiper_doubles.sam" O="STO15_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_43-75" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_hybPiper_doubles.sam" O="STO17_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_33-10" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_hybPiper_doubles.sam" O="STO18_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_560-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_hybPiper_doubles.sam" O="STO19_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_528-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_hybPiper_doubles.sam" O="STO20_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_559-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_hybPiper_doubles.sam" O="STO21_hybPiper_IDed.sam" RGPU="STO" RGSM="STO_33-272" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_hybPiper_doubles.sam" O="SUM01_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_41-103" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_hybPiper_doubles.sam" O="SUM02_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_41-104" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_hybPiper_doubles.sam" O="SUM03_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RA645" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_hybPiper_doubles.sam" O="SUM04_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP31" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_hybPiper_doubles.sam" O="SUM05_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_53-255" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_hybPiper_doubles.sam" O="SUM06_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_53-252" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_hybPiper_doubles.sam" O="SUM07_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_53-253" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_hybPiper_doubles.sam" O="SUM08_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_53-256" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_hybPiper_doubles.sam" O="SUM09_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_hybPiper_doubles.sam" O="SUM10_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_644-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_hybPiper_doubles.sam" O="SUM11_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_hybPiper_doubles.sam" O="SUM12_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_hybPiper_doubles.sam" O="SUM13_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_hybPiper_doubles.sam" O="SUM15_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_hybPiper_doubles.sam" O="SUM16_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_hybPiper_doubles.sam" O="SUM17_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_624-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_hybPiper_doubles.sam" O="SUM18_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_625-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_hybPiper_doubles.sam" O="SUM19_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_626-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_hybPiper_doubles.sam" O="SUM20_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_643-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_hybPiper_doubles.sam" O="SUM21_hybPiper_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL="illumina" RGLB="SeqCap2012"
#  
#  cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts
#  mv HTAdapter1_transcripts_doubles.sam ELF01_transcripts_doubles.sam
#  mv HTAdapter2_transcripts_doubles.sam ELF02_transcripts_doubles.sam
#  mv HTAdapter3_transcripts_doubles.sam ELF03_transcripts_doubles.sam
#  mv HTAdapter4_transcripts_doubles.sam ELF04_transcripts_doubles.sam
#  mv HTAdapter5_transcripts_doubles.sam ELF05_transcripts_doubles.sam
#  mv HTAdapter6_transcripts_doubles.sam ELF06_transcripts_doubles.sam
#  mv HTAdapter7_transcripts_doubles.sam ELF07_transcripts_doubles.sam
#  mv HTAdapter8_transcripts_doubles.sam ELF08_transcripts_doubles.sam
#  mv HTAdapter9_transcripts_doubles.sam ELF09_transcripts_doubles.sam
#  mv HTAdapter10_transcripts_doubles.sam ELF10_transcripts_doubles.sam
#  mv HTAdapter11_transcripts_doubles.sam ELF11_transcripts_doubles.sam
#  mv HTAdapter12_transcripts_doubles.sam ELF12_transcripts_doubles.sam
#  mv HTAdapter13_transcripts_doubles.sam ELF13_transcripts_doubles.sam
#  mv HTAdapter14_transcripts_doubles.sam ELF14_transcripts_doubles.sam
#  mv HTAdapter15_transcripts_doubles.sam ELF15_transcripts_doubles.sam
#  mv HTAdapter16_transcripts_doubles.sam ELF16_transcripts_doubles.sam
#  mv HTAdapter17_transcripts_doubles.sam MAH01_transcripts_doubles.sam
#  mv HTAdapter18_transcripts_doubles.sam MAH02_transcripts_doubles.sam
#  mv HTAdapter19_transcripts_doubles.sam MAH03_transcripts_doubles.sam
#  mv HTAdapter20_transcripts_doubles.sam MAH04_transcripts_doubles.sam
#  mv HTAdapter21_transcripts_doubles.sam MAH05_transcripts_doubles.sam
#  mv HTAdapter22_transcripts_doubles.sam MAH06_transcripts_doubles.sam
#  mv HTAdapter23_transcripts_doubles.sam MAH07_transcripts_doubles.sam
#  mv HTAdapter24_transcripts_doubles.sam MAH08_transcripts_doubles.sam
#  mv HTAdapter25_transcripts_doubles.sam MAH09_transcripts_doubles.sam
#  mv HTAdapter26_transcripts_doubles.sam MAH10_transcripts_doubles.sam
#  mv HTAdapter27_transcripts_doubles.sam MAH11_transcripts_doubles.sam
#  mv HTAdapter28_transcripts_doubles.sam MAH12_transcripts_doubles.sam
#  mv HTAdapter29_transcripts_doubles.sam MAH13_transcripts_doubles.sam
#  mv HTAdapter30_transcripts_doubles.sam MAH14_transcripts_doubles.sam
#  mv HTAdapter31_transcripts_doubles.sam MAH15_transcripts_doubles.sam
#  mv HTAdapter32_transcripts_doubles.sam MAR01_transcripts_doubles.sam
#  mv HTAdapter33_transcripts_doubles.sam MAR02_transcripts_doubles.sam
#  mv HTAdapter34_transcripts_doubles.sam MAR03_transcripts_doubles.sam
#  mv HTAdapter35_transcripts_doubles.sam MAR04_transcripts_doubles.sam
#  mv HTAdapter36_transcripts_doubles.sam MAR05_transcripts_doubles.sam
#  mv HTAdapter37_transcripts_doubles.sam MAR06_transcripts_doubles.sam
#  mv HTAdapter38_transcripts_doubles.sam MAR07_transcripts_doubles.sam
#  mv HTAdapter39_transcripts_doubles.sam MAR08_transcripts_doubles.sam
#  mv HTAdapter40_transcripts_doubles.sam MAR09_transcripts_doubles.sam
#  mv HTAdapter41_transcripts_doubles.sam MAR10_transcripts_doubles.sam
#  mv HTAdapter42_transcripts_doubles.sam PAP01_transcripts_doubles.sam
#  mv HTAdapter43_transcripts_doubles.sam PAP02_transcripts_doubles.sam
#  mv HTAdapter44_transcripts_doubles.sam PAP03_transcripts_doubles.sam
#  mv HTAdapter45_transcripts_doubles.sam PAP04_transcripts_doubles.sam
#  mv HTAdapter46_transcripts_doubles.sam PAP05_transcripts_doubles.sam
#  mv HTAdapter47_transcripts_doubles.sam PAP06_transcripts_doubles.sam
#  mv HTAdapter48_transcripts_doubles.sam PAP07_transcripts_doubles.sam
#  mv HTAdapter49_transcripts_doubles.sam PAP08_transcripts_doubles.sam
#  mv HTAdapter50_transcripts_doubles.sam PAP09_transcripts_doubles.sam
#  mv HTAdapter51_transcripts_doubles.sam PAP10_transcripts_doubles.sam
#  mv HTAdapter52_transcripts_doubles.sam PAP11_transcripts_doubles.sam
#  mv HTAdapter53_transcripts_doubles.sam PAP12_transcripts_doubles.sam
#  mv HTAdapter54_transcripts_doubles.sam PAP13_transcripts_doubles.sam
#  mv HTAdapter55_transcripts_doubles.sam PAP14_transcripts_doubles.sam
#  mv HTAdapter56_transcripts_doubles.sam PAP15_transcripts_doubles.sam
#  mv HTAdapter57_transcripts_doubles.sam STO01_transcripts_doubles.sam
#  mv HTAdapter58_transcripts_doubles.sam STO02_transcripts_doubles.sam
#  mv HTAdapter59_transcripts_doubles.sam STO03_transcripts_doubles.sam
#  mv HTAdapter60_transcripts_doubles.sam STO04_transcripts_doubles.sam
#  mv HTAdapter61_transcripts_doubles.sam STO05_transcripts_doubles.sam
#  mv HTAdapter62_transcripts_doubles.sam STO06_transcripts_doubles.sam
#  mv HTAdapter63_transcripts_doubles.sam STO07_transcripts_doubles.sam
#  mv HTAdapter64_transcripts_doubles.sam STO08_transcripts_doubles.sam
#  mv HTAdapter65_transcripts_doubles.sam STO09_transcripts_doubles.sam
#  mv HTAdapter66_transcripts_doubles.sam STO10_transcripts_doubles.sam
#  mv HTAdapter67_transcripts_doubles.sam STO11_transcripts_doubles.sam
#  mv HTAdapter68_transcripts_doubles.sam STO12_transcripts_doubles.sam
#  mv HTAdapter69_transcripts_doubles.sam STO13_transcripts_doubles.sam
#  mv HTAdapter70_transcripts_doubles.sam STO14_transcripts_doubles.sam
#  mv HTAdapter71_transcripts_doubles.sam STO15_transcripts_doubles.sam
#  mv HTAdapter72_transcripts_doubles.sam STO17_transcripts_doubles.sam
#  mv HTAdapter73_transcripts_doubles.sam STO18_transcripts_doubles.sam
#  mv HTAdapter74_transcripts_doubles.sam STO19_transcripts_doubles.sam
#  mv HTAdapter75_transcripts_doubles.sam STO20_transcripts_doubles.sam
#  mv HTAdapter76_transcripts_doubles.sam STO21_transcripts_doubles.sam
#  mv HTAdapter77_transcripts_doubles.sam SUM01_transcripts_doubles.sam
#  mv HTAdapter78_transcripts_doubles.sam SUM02_transcripts_doubles.sam
#  mv HTAdapter79_transcripts_doubles.sam SUM03_transcripts_doubles.sam
#  mv HTAdapter80_transcripts_doubles.sam SUM04_transcripts_doubles.sam
#  mv HTAdapter81_transcripts_doubles.sam SUM05_transcripts_doubles.sam
#  mv HTAdapter82_transcripts_doubles.sam SUM06_transcripts_doubles.sam
#  mv HTAdapter83_transcripts_doubles.sam SUM07_transcripts_doubles.sam
#  mv HTAdapter84_transcripts_doubles.sam SUM08_transcripts_doubles.sam
#  mv HTAdapter85_transcripts_doubles.sam SUM09_transcripts_doubles.sam
#  mv HTAdapter86_transcripts_doubles.sam SUM10_transcripts_doubles.sam
#  mv HTAdapter87_transcripts_doubles.sam SUM11_transcripts_doubles.sam
#  mv HTAdapter88_transcripts_doubles.sam SUM12_transcripts_doubles.sam
#  mv HTAdapter89_transcripts_doubles.sam SUM13_transcripts_doubles.sam
#  mv HTAdapter90_transcripts_doubles.sam SUM15_transcripts_doubles.sam
#  mv HTAdapter91_transcripts_doubles.sam SUM16_transcripts_doubles.sam
#  mv HTAdapter92_transcripts_doubles.sam SUM17_transcripts_doubles.sam
#  mv HTAdapter93_transcripts_doubles.sam SUM18_transcripts_doubles.sam
#  mv HTAdapter94_transcripts_doubles.sam SUM19_transcripts_doubles.sam
#  mv HTAdapter95_transcripts_doubles.sam SUM20_transcripts_doubles.sam
#  mv HTAdapter96_transcripts_doubles.sam SUM21_transcripts_doubles.sam
#  
#  ## -- ADD READGROUPS -- ##
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_transcripts_doubles.sam" O="ELF01_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_54-38" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_transcripts_doubles.sam" O="ELF02_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RA607" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_transcripts_doubles.sam" O="ELF03_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_transcripts_doubles.sam" O="ELF04_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_transcripts_doubles.sam" O="ELF05_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_558-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_transcripts_doubles.sam" O="ELF06_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_565-13" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_transcripts_doubles.sam" O="ELF07_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_5651-15" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_transcripts_doubles.sam" O="ELF08_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_568-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_transcripts_doubles.sam" O="ELF09_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_546-03" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_transcripts_doubles.sam" O="ELF10_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_547-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_transcripts_doubles.sam" O="ELF11_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_6050-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_transcripts_doubles.sam" O="ELF12_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_622-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_transcripts_doubles.sam" O="ELF13_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_636-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_transcripts_doubles.sam" O="ELF14_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_transcripts_doubles.sam" O="ELF15_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_transcripts_doubles.sam" O="ELF16_transcripts_IDed.sam" RGPU="ELF" RGSM="ELF_RP54" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_transcripts_doubles.sam" O="MAH01_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_35-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_transcripts_doubles.sam" O="MAH02_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_520-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_transcripts_doubles.sam" O="MAH03_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_539-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_transcripts_doubles.sam" O="MAH04_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_540-07" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_transcripts_doubles.sam" O="MAH05_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_541-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_transcripts_doubles.sam" O="MAH06_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_542-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_transcripts_doubles.sam" O="MAH07_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_543-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_transcripts_doubles.sam" O="MAH08_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RA620" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_transcripts_doubles.sam" O="MAH09_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RA641" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_transcripts_doubles.sam" O="MAH10_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_transcripts_doubles.sam" O="MAH11_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_transcripts_doubles.sam" O="MAH12_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_transcripts_doubles.sam" O="MAH13_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_transcripts_doubles.sam" O="MAH14_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_transcripts_doubles.sam" O="MAH15_transcripts_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_transcripts_doubles.sam" O="MAR01_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RP6771" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_transcripts_doubles.sam" O="MAR02_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RP686" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_transcripts_doubles.sam" O="MAR03_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_368- 2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_transcripts_doubles.sam" O="MAR04_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_transcripts_doubles.sam" O="MAR05_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA379" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_transcripts_doubles.sam" O="MAR06_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_transcripts_doubles.sam" O="MAR07_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_transcripts_doubles.sam" O="MAR08_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA605" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_transcripts_doubles.sam" O="MAR09_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA630" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_transcripts_doubles.sam" O="MAR10_transcripts_IDed.sam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_transcripts_doubles.sam" O="PAP01_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_561-4" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_transcripts_doubles.sam" O="PAP02_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_562-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_transcripts_doubles.sam" O="PAP03_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_564-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_transcripts_doubles.sam" O="PAP04_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_569-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_transcripts_doubles.sam" O="PAP05_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_5010-05" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_transcripts_doubles.sam" O="PAP06_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_5021-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_transcripts_doubles.sam" O="PAP07_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_505-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_transcripts_doubles.sam" O="PAP08_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_516-02" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_transcripts_doubles.sam" O="PAP09_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_5241-01" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_transcripts_doubles.sam" O="PAP10_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_28-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_transcripts_doubles.sam" O="PAP11_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_RA434" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_transcripts_doubles.sam" O="PAP12_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_RA647" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_transcripts_doubles.sam" O="PAP13_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_RA648" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_transcripts_doubles.sam" O="PAP14_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_RA649" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_transcripts_doubles.sam" O="PAP15_transcripts_IDed.sam" RGPU="PAP" RGSM="PAP_RP16" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_transcripts_doubles.sam" O="STO01_transcripts_IDed.sam" RGPU="STO" RGSM="STO_RP697" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_transcripts_doubles.sam" O="STO02_transcripts_IDed.sam" RGPU="STO" RGSM="STO_RP698" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_transcripts_doubles.sam" O="STO03_transcripts_IDed.sam" RGPU="STO" RGSM="STO_RA530" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_transcripts_doubles.sam" O="STO04_transcripts_IDed.sam" RGPU="STO" RGSM="STO_RA549" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_transcripts_doubles.sam" O="STO05_transcripts_IDed.sam" RGPU="STO" RGSM="STO_RP22" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_transcripts_doubles.sam" O="STO06_transcripts_IDed.sam" RGPU="STO" RGSM="STO_271" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_transcripts_doubles.sam" O="STO07_transcripts_IDed.sam" RGPU="STO" RGSM="STO_275" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_transcripts_doubles.sam" O="STO08_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-113" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_transcripts_doubles.sam" O="STO09_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-115" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_transcripts_doubles.sam" O="STO10_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-125" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_transcripts_doubles.sam" O="STO11_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-76" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_transcripts_doubles.sam" O="STO12_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-77" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_transcripts_doubles.sam" O="STO13_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-78" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_transcripts_doubles.sam" O="STO14_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-74" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_transcripts_doubles.sam" O="STO15_transcripts_IDed.sam" RGPU="STO" RGSM="STO_43-75" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_transcripts_doubles.sam" O="STO17_transcripts_IDed.sam" RGPU="STO" RGSM="STO_33-10" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_transcripts_doubles.sam" O="STO18_transcripts_IDed.sam" RGPU="STO" RGSM="STO_560-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_transcripts_doubles.sam" O="STO19_transcripts_IDed.sam" RGPU="STO" RGSM="STO_528-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_transcripts_doubles.sam" O="STO20_transcripts_IDed.sam" RGPU="STO" RGSM="STO_559-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_transcripts_doubles.sam" O="STO21_transcripts_IDed.sam" RGPU="STO" RGSM="STO_33-272" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_transcripts_doubles.sam" O="SUM01_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_41-103" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_transcripts_doubles.sam" O="SUM02_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_41-104" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_transcripts_doubles.sam" O="SUM03_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RA645" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_transcripts_doubles.sam" O="SUM04_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP31" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_transcripts_doubles.sam" O="SUM05_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_53-255" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_transcripts_doubles.sam" O="SUM06_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_53-252" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_transcripts_doubles.sam" O="SUM07_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_53-253" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_transcripts_doubles.sam" O="SUM08_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_53-256" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_transcripts_doubles.sam" O="SUM09_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_transcripts_doubles.sam" O="SUM10_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_644-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_transcripts_doubles.sam" O="SUM11_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_transcripts_doubles.sam" O="SUM12_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_transcripts_doubles.sam" O="SUM13_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_transcripts_doubles.sam" O="SUM15_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_transcripts_doubles.sam" O="SUM16_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_transcripts_doubles.sam" O="SUM17_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_624-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_transcripts_doubles.sam" O="SUM18_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_625-1" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_transcripts_doubles.sam" O="SUM19_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_626-3" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_transcripts_doubles.sam" O="SUM20_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_643-2" RGPL="illumina" RGLB="SeqCap2012"
#  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_transcripts_doubles.sam" O="SUM21_transcripts_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL="illumina" RGLB="SeqCap2012"
#  echo "RLK_report: readGrouping complete"
# 
# # RNA Seq
# cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629651_doubles.sam" O="SRR629651_IDed.sam" RGPU="SUM" RGSM="SAMN01823448" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629599_doubles.sam" O="SRR629599_IDed.sam" RGPU="SUM" RGSM="SAMN01823447" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629652_doubles.sam" O="SRR629652_IDed.sam" RGPU="SUM" RGSM="SAMN01823449" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629653_doubles.sam" O="SRR629653_IDed.sam" RGPU="SUM" RGSM="SAMN01823450" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629654_doubles.sam" O="SRR629654_IDed.sam" RGPU="SUM" RGSM="SAMN01823451" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629655_doubles.sam" O="SRR629655_IDed.sam" RGPU="SUM" RGSM="SAMN01823452" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629656_doubles.sam" O="SRR629656_IDed.sam" RGPU="SUM" RGSM="SAMN01823453" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629657_doubles.sam" O="SRR629657_IDed.sam" RGPU="SUM" RGSM="SAMN01823454" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629658_doubles.sam" O="SRR629658_IDed.sam" RGPU="SUM" RGSM="SAMN01823455" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629659_doubles.sam" O="SRR629659_IDed.sam" RGPU="SUM" RGSM="SAMN01823456" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629661_doubles.sam" O="SRR629661_IDed.sam" RGPU="SUM" RGSM="SAMN01823457" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629660_doubles.sam" O="SRR629660_IDed.sam" RGPU="SUM" RGSM="SAMN01823458" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629662_doubles.sam" O="SRR629662_IDed.sam" RGPU="SUM" RGSM="SAMN01823459" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629664_doubles.sam" O="SRR629664_IDed.sam" RGPU="SUM" RGSM="SAMN01823460" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629663_doubles.sam" O="SRR629663_IDed.sam" RGPU="SUM" RGSM="SAMN01823461" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629665_doubles.sam" O="SRR629665_IDed.sam" RGPU="SUM" RGSM="SAMN01823462" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629666_doubles.sam" O="SRR629666_IDed.sam" RGPU="SUM" RGSM="SAMN01823463" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629667_doubles.sam" O="SRR629667_IDed.sam" RGPU="SUM" RGSM="SAMN01823464" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629668_doubles.sam" O="SRR629668_IDed.sam" RGPU="SUM" RGSM="SAMN01823465" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629669_doubles.sam" O="SRR629669_IDed.sam" RGPU="SUM" RGSM="SAMN01823466" RGPL="illumina" RGLB="RNASeq2012" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497737_doubles.sam" O="SRR497737_IDed.sam" RGPU="SUM" RGSM="SAMN00996377" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497739_doubles.sam" O="SRR497739_IDed.sam" RGPU="SUM" RGSM="SAMN00996379" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497740_doubles.sam" O="SRR497740_IDed.sam" RGPU="SUM" RGSM="SAMN00996380" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497741_doubles.sam" O="SRR497741_IDed.sam" RGPU="SUM" RGSM="SAMN00996381" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497742_doubles.sam" O="SRR497742_IDed.sam" RGPU="SUM" RGSM="SAMN00996382" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497743_doubles.sam" O="SRR497743_IDed.sam" RGPU="SUM" RGSM="SAMN00996383" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497744_doubles.sam" O="SRR497744_IDed.sam" RGPU="SUM" RGSM="SAMN00996384" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497745_doubles.sam" O="SRR497745_IDed.sam" RGPU="SUM" RGSM="SAMN00996385" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497746_doubles.sam" O="SRR497746_IDed.sam" RGPU="SUM" RGSM="SAMN00996386" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497747_doubles.sam" O="SRR497747_IDed.sam" RGPU="SUM" RGSM="SAMN00996387" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497748_doubles.sam" O="SRR497748_IDed.sam" RGPU="SUM" RGSM="SAMN00996388" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497749_doubles.sam" O="SRR497749_IDed.sam" RGPU="SUM" RGSM="SAMN00996389" RGPL="illumina" RGLB="RNASeq2008" 
# java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497738_doubles.sam" O="SRR497738_IDed.sam" RGPU="SUM" RGSM="SAMN00996378" RGPL="illumina" RGLB="RNASeq2008" 

# -- PREPARE FOR SNP CALLING -- ##

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ HYB PIPER BLOCK ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# HybPiper
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper
# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# make result directories
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperSNPTables
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperSequenceTables
# copy reference for SNP calling
cp /home/rlk0015/SeqCap/code/References/HybPiper_Transcripts_Round1.fasta /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa
# index ref
samtools faidx /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa
java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
    R= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    O= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.dict
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper
while read i;
do
# convert .sam to .bam & sort ###
samtools view -@ 2 -bS /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper/"$i"_sorted.bam
### remove duplicates ###
java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
    INPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/HybPiper/"$i"_sorted.bam \
    OUTPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_dupsRemoved.bam \
    METRICS_FILE=DuplicationMetrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORTED=TRUE \
    REMOVE_DUPLICATES=TRUE
#index sorted bam
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_dupsRemoved.bam

## Calculate Stats ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory_ReducedTranscriptome/mappedReadsDNA/HybPiper
# make stats folder
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperStats
# tally mapped reads & calcuate the stats
samtools idxstats /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperStats/"$i"_counts.txt
samtools flagstat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperStats/"$i"_stats.txt
samtools depth /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperStats/"$i"_depth.txt
done<samList

cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperStats
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList

# Generate file with all depth information
cat *_depth.txt > HybPiperContigs_depth.txt
for f in HybPiperContigs_depth.txt
do
sed -i "s/$/\tHybPiperContigs/" $f; done

# Create file with avg depth per exon using Randy's python script
#!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py HybPiperContigs_depth.txt HybPiperContigs_avgDepth.txt

#  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
#  # make results directory & move results
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/stats
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/counts
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/avgDepth
#  
#  # cp results to respective directories
#  cp *stats.txt /home/rlk0015/SeqCap/Dec2018/stats
#  cp *counts.txt /home/rlk0015/SeqCap/Dec2018/counts
#  cp *depth.txt /home/rlk0015/SeqCap/Dec2018/avgDepth
#  *****************************************************************************************************


### move to GATK directory ###
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/
### merge .bam files ###
samtools merge /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/dupsRemoved.bam /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/*_dupsRemoved.bam
### index the merged .bam ###
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/dupsRemoved.bam
# call indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/dupsRemoved.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsCalled.intervals
# realign indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/dupsRemoved.bam \
    -targetIntervals /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsCalled.intervals \
    -LOD 3.0 \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsRealigned.bam
## -- CALL SNPS -- ##
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsRealigned.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/calledSNPs.vcf \
    -gt_mode DISCOVERY \
    -ploidy 2 \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# annotate variants
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsRealigned.bam \
    -G StandardAnnotation \
    -V:variant,VCF /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/calledSNPs.vcf \
    -XA SnpEff \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/annotatedSNPs.vcf \
    -rf BadCigar
# annotate indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsRealigned.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/annotatedIndels.vcf \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# mask indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/calledSNPs.vcf \
    --mask /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/annotatedIndels.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/SNPsMaskedIndels.vcf/ \
    -rf BadCigar
# restrict to high-quality variant calls
cat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/qualitySNPs.vcf
# read-backed phasing
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/indelsRealigned.bam \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/qualitySNPs.vcf \
    -L /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/qualitySNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/phasedSNPs.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar

## Make Sample List from VCF ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK
bcftools query -l phasedSNPs.vcf > VcfSampleList

while read i;
do
# VCF for each sample
java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/phasedSNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_phasedSNPs.vcf \
    -sn "$i" \
    -rf BadCigar
# make SNPs table
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/"$i"_phasedSNPs.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperSNPTables/"$i"_tableSNPs.txt \
    -rf BadCigar
# Add phased SNPs to reference and filter
python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperGATK/HybPiperContigs.fa \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperSNPTables/"$i"_tableSNPs.txt \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/HybPiperSequenceTables/"$i"_tableSequences.txt \
    1
done<VcfSampleList

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ TRANSCRIPTS BLOCK ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts
# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# make result directories
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsSNPTables
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsSequenceTables
# copy reference for SNP calling
cp /home/rlk0015/SeqCap/code/References/Transcripts.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa
# index ref
samtools faidx /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa
java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
    R= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    O= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.dict
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts
while read i;
do
# convert .sam to .bam & sort ###
samtools view -@ 2 -bS /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"_sorted.bam
### remove duplicates ###
java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
    INPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"_sorted.bam \
    OUTPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam \
    METRICS_FILE=DuplicationMetrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORTED=TRUE \
    REMOVE_DUPLICATES=TRUE
#index sorted bam
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam

## Calculate Stats ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory_ReducedTranscriptome/mappedReadsDNA/Transcripts
# make stats folder
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsStats
# tally mapped reads & calcuate the stats
samtools idxstats /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsStats/"$i"_counts.txt
samtools flagstat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsStats/"$i"_stats.txt
samtools depth /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsStats/"$i"_depth.txt
done<samList

cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsStats
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList

# Generate file with all depth information
cat *_depth.txt > Transcripts_depth.txt
for f in Transcripts_depth.txt
do
sed -i "s/$/\tTranscripts/" $f; done

# Create file with avg depth per exon using Randy's python script
#!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py Transcripts_depth.txt Transcripts_avgDepth.txt

#  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
#  # make results directory & move results
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/stats
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/counts
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/avgDepth
#  
#  # cp results to respective directories
#  cp *stats.txt /home/rlk0015/SeqCap/Dec2018/stats
#  cp *counts.txt /home/rlk0015/SeqCap/Dec2018/counts
#  cp *depth.txt /home/rlk0015/SeqCap/Dec2018/avgDepth
#  *****************************************************************************************************


### move to GATK directory ###
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/
### merge .bam files ###
samtools merge /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/*_dupsRemoved.bam
### index the merged .bam ###
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam
# call indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsCalled.intervals
# realign indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam \
    -targetIntervals /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsCalled.intervals \
    -LOD 3.0 \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam
## -- CALL SNPS -- ##
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/calledSNPs.vcf \
    -gt_mode DISCOVERY \
    -ploidy 2 \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# annotate variants
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
    -G StandardAnnotation \
    -V:variant,VCF /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/calledSNPs.vcf \
    -XA SnpEff \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/annotatedSNPs.vcf \
    -rf BadCigar
# annotate indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/annotatedIndels.vcf \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# mask indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/calledSNPs.vcf \
    --mask /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/annotatedIndels.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/SNPsMaskedIndels.vcf/ \
    -rf BadCigar
# restrict to high-quality variant calls
cat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf
# read-backed phasing
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf \
    -L /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/phasedSNPs.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar

## Make Sample List from VCF ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK
bcftools query -l phasedSNPs.vcf > VcfSampleList

while read i;
do
# VCF for each sample
java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/phasedSNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_phasedSNPs.vcf \
    -sn "$i" \
    -rf BadCigar
# make SNPs table
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/"$i"_phasedSNPs.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsSNPTables/"$i"_tableSNPs.txt \
    -rf BadCigar
# Add phased SNPs to reference and filter
python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsSNPTables/"$i"_tableSNPs.txt \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/TranscriptsSequenceTables/"$i"_tableSequences.txt \
    1
done<VcfSampleList

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ RNA SEQ BLOCK ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# make result directories
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNASNPTables
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNASequenceTables
# copy reference for SNP calling
cp /home/rlk0015/SeqCap/code/References/TranscriptomePlusTranscripts.fa /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa
# index ref
samtools faidx /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa
java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
    R= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    O= /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.dict
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA
while read i;
do
# convert .sam to .bam & sort ###
samtools view -@ 2 -bS /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam
### remove duplicates ###
java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
    INPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam \
    OUTPUT=/scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam \
    METRICS_FILE=DuplicationMetrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORTED=TRUE \
    REMOVE_DUPLICATES=TRUE
#index sorted bam
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam

## Calculate Stats ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory_ReducedTranscriptome/mappedReadsRNA
# make stats folder
mkdir -p /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAStats
# tally mapped reads & calcuate the stats
samtools idxstats /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAStats/"$i"_counts.txt
samtools flagstat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAStats/"$i"_stats.txt
samtools depth /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAStats/"$i"_depth.txt
done<samList

cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAStats
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList

# Generate file with all depth information
cat *_depth.txt > RNA_depth.txt
for f in RNA_depth.txt
do
sed -i "s/$/\tRNA/" $f; done

# Create file with avg depth per exon using Randy's python script
#!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py RNA_depth.txt RNA_avgDepth.txt

#  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
#  # make results directory & move results
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/stats
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/counts
#  mkdir -p /home/rlk0015/SeqCap/Dec2018/avgDepth
#  
#  # cp results to respective directories
#  cp *stats.txt /home/rlk0015/SeqCap/Dec2018/stats
#  cp *counts.txt /home/rlk0015/SeqCap/Dec2018/counts
#  cp *depth.txt /home/rlk0015/SeqCap/Dec2018/avgDepth
#  *****************************************************************************************************


### move to GATK directory ###
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/
### merge .bam files ###
samtools merge /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/dupsRemoved.bam /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/*_dupsRemoved.bam
### index the merged .bam ###
samtools index /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/dupsRemoved.bam
# call indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/dupsRemoved.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsCalled.intervals
# realign indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/dupsRemoved.bam \
    -targetIntervals /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsCalled.intervals \
    -LOD 3.0 \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsRealigned.bam
## -- CALL SNPS -- ##
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsRealigned.bam \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/calledSNPs.vcf \
    -gt_mode DISCOVERY \
    -ploidy 2 \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# annotate variants
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsRealigned.bam \
    -G StandardAnnotation \
    -V:variant,VCF /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/calledSNPs.vcf \
    -XA SnpEff \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/annotatedSNPs.vcf \
    -rf BadCigar
# annotate indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsRealigned.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/annotatedIndels.vcf \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# mask indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/calledSNPs.vcf \
    --mask /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/annotatedIndels.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/SNPsMaskedIndels.vcf/ \
    -rf BadCigar
# restrict to high-quality variant calls
cat /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/qualitySNPs.vcf
# read-backed phasing
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -I /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/indelsRealigned.bam \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/qualitySNPs.vcf \
    -L /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/qualitySNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/phasedSNPs.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar

## Make Sample List from VCF ##
cd /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK
bcftools query -l phasedSNPs.vcf > VcfSampleList

while read i;
do
# VCF for each sample
java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    --variant /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/phasedSNPs.vcf \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_phasedSNPs.vcf \
    -sn "$i" \
    -rf BadCigar
# make SNPs table
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    -V /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/"$i"_phasedSNPs.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNASNPTables/"$i"_tableSNPs.txt \
    -rf BadCigar
# Add phased SNPs to reference and filter
python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNAGATK/Transcriptome.fa \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNASNPTables/"$i"_tableSNPs.txt \
    /scratch/rlk0015/Telag/Dec2018/WorkingDirectory/RNASequenceTables/"$i"_tableSequences.txt \
    1
done<VcfSampleList
