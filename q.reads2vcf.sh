#!/bin/sh

#Give job a name
#PBS -N FullScript_May2020

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=1:ppn=5,walltime=40:00:00 
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
module load bcftools/1.3.2
module load hybpiper/1
module load xz/5.2.2
module load htslib/1.3.3
module load python/3.6.4
module load gatk/4.1.2.0

#create working directory ###
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory
# DNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/assembledReads
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsDNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA
  # -- make result directories
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/StatsDNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/SNPTablesDNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/SequenceTablesDNA
# RNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsRNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA
# Utilities
mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/References
echo "RLK_report: directory created: /scratch/rlk0015/Telag/May2020/WorkingDirectory with rawReads and cleanReads sub directories"

# ------ For optimization purposes we'll just work with a couple of individuals, uncomment when performing complete script ----------#
# ### copy raw data to scratch ###
# cd /home/shared/tss0019_lab/SeqCap_GarterSnake2012/
# for i in {1..96}
# do
#   cp Sample_HTAdapter"$i"/*.fastq.gz /scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsDNA
# done
# -----------------------------------------------------------------------------------------------------------------------------------#
#+  +++++++++++++ Done May 2020 +++++++++++++++++ #
#+  ### perform initial quality check ###
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsDNA
#+  ls *.fastq.gz | time parallel -j+0 --eta 'fastqc {}'
#+  multiqc .
#+  
#+  ### copy over adapter file ###
#+  cp /home/rlk0015/SeqCap/code/References/adapters.fa /scratch/rlk0015/Telag/May2020/WorkingDirectory/References
#+  
#+  ### paired-end trimming ###
#+  ls | grep "fastq.gz" | cut -d "_" -f 1,2 | sort | uniq > PE_TrimmList
#+  while read i
#+  do
#+  java -jar /tools/trimmomatic-0.36/trimmomatic-0.36.jar \
#+  	PE \
#+  	-threads 6 \
#+  	-phred33 \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsDNA/"$i"*_R1*.fastq.gz \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/rawReadsDNA/"$i"*_R2*.fastq.gz \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_R1_paired.fastq.gz \ 
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_R1_unpaired.fastq.gz \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_R2_paired.fastq.gz \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_R2_unpaired.fastq.gz \
#+  	ILLUMINACLIP:adapters.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:30 MINLEN:36
#+  done<PE_TrimmList
#+  ### perform second quality check ###
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/
#+  ls *paired.fastq.gz | grep "paired.fastq.gz" | time parallel -j+0 --eta 'fastqc {}'
#+  multiqc .
#+  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---------------------------
# Copy clean RNA reads over
# ---------------------------
# ------ For optimization purposes we'll just work with a couple of individuals, uncomment when performing complete script ----------#
# cd /scratch/GarterSnake/RNAseq_2012Data/CleanedData
# ls SRR*.fastq | parallel -j+0 --eta 'cp {} /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA'
# mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/Report/RNAseq_2012copied
# cd /scratch/GarterSnake/RNAseq_2008Data/CleanedData
# ls SRR*.fastq | parallel -j+0 --eta 'cp {} /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA'
# mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/Report/RNAseq_2008copied
# -----------------------------------------------------------------------------------------------------------------------------------#
#### ********************************************* ####

# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººTonia said to not worry about HybPiper for nowººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººBut this is a working script if we change mindsºººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# # -------- HybPiper Assembly -------- ###
# cp /home/rlk0015/SeqCap/code/References/Transcripts.fa /scratch/rlk0015/Telag/May2020/WorkingDirectory/assembledReads
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA
# ls | grep "fastq" | cut -d "_" -f 1 | sort | uniq > cleanReadsList
# while read i
# do
# reads_first.py -b /scratch/rlk0015/Telag/May2020/WorkingDirectory/assembledReads/Transcripts.fa -r /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"*_R*_paired.fastq --bwa --cpu 20 --prefix /scratch/rlk0015/Telag/May2020/WorkingDirectory/assembledReads/HybPiperAssembly
# done<cleanReadsList
# # Create test_seq_lengths.txt file (for data visualization)
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/assembledReads
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
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #

# -------------- Already performed- see details below --------------
# # Tonia's criteria for trimming RNA-Seq reads:
# ########## 2008 Data Trimmomatic #############
# java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar SE -phred33 "$i"_1.fastq  /scratch/GarterSnake/RNAseq_2008Data/CleanedData/"$i"_cleaned.fastq  LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:36 2012 samples are PE 100 bp reads.
# ########## 2012 Trimmomatic #############
# java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar PE  -phred33 "$i"_1.fastq "$i"_2.fastq "$i"_1_paired.fastq "$i"_1_unpaired.fastq "$i"_2_paired.fastq "$i"_2_unpaired.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
# ------------------------------------------------------------------

#+  +++++++++++++ Done May 2020 +++++++++++++++++ #
#+  # Quality check on RNA-Seq data
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA
#+  ls *.fastq | grep "fastq" | time parallel -j+0 --eta 'fastqc {}'
#+  multiqc . 
### ***  -------------------  MAP DNA TO REFERENCE  ------------------ *** ###
### ••••••••••••••••••••••• This section is complete ••••••••••••••••••••• ###
# ### copy the reference (T. elegans genome) to References directory ###
# cp /home/rlk0015/SeqCap/code/References/T_elegans_genome/latest_assembly_versions/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.fna.gz /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fna.gz
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/References
# gunzip TelagGenome.fna.gz
# echo "RLK_report: reference copy complete"
### •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• ###
#+  +++++++++++++ Done May 2020 +++++++++++++++++ #
#+  # index T. elegans genome
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/References
#+  bwa index -p TelagGenomeBWA -a is TelagGenome.fna
#+  echo "RLK_report: reference index complete"
#+  # create list with each paired individual
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA
#+  ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/pairedMapList
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA
#+  # while loop through the names in pairedMapList
#+  while read i
#+  do
#+  ### map to T. elegans genome ###
#+  bwa mem -t 4 
#+  -M /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenomeBWA 
#+  scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_*R1_paired.fastq 
#+  /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsDNA/"$i"_*R2_paired.fastq > /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/"$i"_mapped.sam
#+  done</scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/pairedMapList
#+  
#+  ### ***  -------------------  MAP RNA TO REFERENCE  ------------------ *** ###
#+  # create list with each paired individual
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA
#+  ls | grep "_cleaned.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/singleEndMapList
#+  ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/pairedEndMapList
#+  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA
#+  # while loop through the names in pairedEndMapList
#+  while read i
#+  do
#+  ### PAIRED-END MAPPING ###
#+  bwa mem -t 4 \ 
#+  	-M  /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenomeBWA \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/"$i"*_1_paired.fastq \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/"$i"*_2_paired.fastq > \
#+  	/scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA/"$i"_mapped.sam
#+  done</scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/pairedEndMapList
#+  echo "RLK_report: paired map complete"
#+  ### SINGLE-END MAPPING ###
#+  while read i
#+  do
#+  bwa mem -t 4
#+  -M /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenomeBWA
#+  /scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/"$i"*.fastq > /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA/"$i"_mapped.sam
#+  done</scratch/rlk0015/Telag/May2020/WorkingDirectory/cleanReadsRNA/singleEndMapList
#+  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ***********************************************************************

### change names to match population sampling ###
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/
mv HTAdapter1_mapped.sam ELF01_mapped.sam
mv HTAdapter2_mapped.sam ELF02_mapped.sam
mv HTAdapter3_mapped.sam ELF03_mapped.sam
mv HTAdapter4_mapped.sam ELF04_mapped.sam
mv HTAdapter5_mapped.sam ELF05_mapped.sam
mv HTAdapter6_mapped.sam ELF06_mapped.sam
mv HTAdapter7_mapped.sam ELF07_mapped.sam
mv HTAdapter8_mapped.sam ELF08_mapped.sam
mv HTAdapter9_mapped.sam ELF09_mapped.sam
mv HTAdapter10_mapped.sam ELF10_mapped.sam
mv HTAdapter11_mapped.sam ELF11_mapped.sam
mv HTAdapter12_mapped.sam ELF12_mapped.sam
mv HTAdapter13_mapped.sam ELF13_mapped.sam
mv HTAdapter14_mapped.sam ELF14_mapped.sam
mv HTAdapter15_mapped.sam ELF15_mapped.sam
mv HTAdapter16_mapped.sam ELF16_mapped.sam
mv HTAdapter17_mapped.sam MAH01_mapped.sam
mv HTAdapter18_mapped.sam MAH02_mapped.sam
mv HTAdapter19_mapped.sam MAH03_mapped.sam
mv HTAdapter20_mapped.sam MAH04_mapped.sam
mv HTAdapter21_mapped.sam MAH05_mapped.sam
mv HTAdapter22_mapped.sam MAH06_mapped.sam
mv HTAdapter23_mapped.sam MAH07_mapped.sam
mv HTAdapter24_mapped.sam MAH08_mapped.sam
mv HTAdapter25_mapped.sam MAH09_mapped.sam
mv HTAdapter26_mapped.sam MAH10_mapped.sam
mv HTAdapter27_mapped.sam MAH11_mapped.sam
mv HTAdapter28_mapped.sam MAH12_mapped.sam
mv HTAdapter29_mapped.sam MAH13_mapped.sam
mv HTAdapter30_mapped.sam MAH14_mapped.sam
mv HTAdapter31_mapped.sam MAH15_mapped.sam
mv HTAdapter32_mapped.sam MAR01_mapped.sam
mv HTAdapter33_mapped.sam MAR02_mapped.sam
mv HTAdapter34_mapped.sam MAR03_mapped.sam
mv HTAdapter35_mapped.sam MAR04_mapped.sam
mv HTAdapter36_mapped.sam MAR05_mapped.sam
mv HTAdapter37_mapped.sam MAR06_mapped.sam
mv HTAdapter38_mapped.sam MAR07_mapped.sam
mv HTAdapter39_mapped.sam MAR08_mapped.sam
mv HTAdapter40_mapped.sam MAR09_mapped.sam
mv HTAdapter41_mapped.sam MAR10_mapped.sam
mv HTAdapter42_mapped.sam PAP01_mapped.sam
mv HTAdapter43_mapped.sam PAP02_mapped.sam
mv HTAdapter44_mapped.sam PAP03_mapped.sam
mv HTAdapter45_mapped.sam PAP04_mapped.sam
mv HTAdapter46_mapped.sam PAP05_mapped.sam
mv HTAdapter47_mapped.sam PAP06_mapped.sam
mv HTAdapter48_mapped.sam PAP07_mapped.sam
mv HTAdapter49_mapped.sam PAP08_mapped.sam
mv HTAdapter50_mapped.sam PAP09_mapped.sam
mv HTAdapter51_mapped.sam PAP10_mapped.sam
mv HTAdapter52_mapped.sam PAP11_mapped.sam
mv HTAdapter53_mapped.sam PAP12_mapped.sam
mv HTAdapter54_mapped.sam PAP13_mapped.sam
mv HTAdapter55_mapped.sam PAP14_mapped.sam
mv HTAdapter56_mapped.sam PAP15_mapped.sam
mv HTAdapter57_mapped.sam STO01_mapped.sam
mv HTAdapter58_mapped.sam STO02_mapped.sam
mv HTAdapter59_mapped.sam STO03_mapped.sam
mv HTAdapter60_mapped.sam STO04_mapped.sam
mv HTAdapter61_mapped.sam STO05_mapped.sam
mv HTAdapter62_mapped.sam STO06_mapped.sam
mv HTAdapter63_mapped.sam STO07_mapped.sam
mv HTAdapter64_mapped.sam STO08_mapped.sam
mv HTAdapter65_mapped.sam STO09_mapped.sam
mv HTAdapter66_mapped.sam STO10_mapped.sam
mv HTAdapter67_mapped.sam STO11_mapped.sam
mv HTAdapter68_mapped.sam STO12_mapped.sam
mv HTAdapter69_mapped.sam STO13_mapped.sam
mv HTAdapter70_mapped.sam STO14_mapped.sam
mv HTAdapter71_mapped.sam STO15_mapped.sam
mv HTAdapter72_mapped.sam STO17_mapped.sam
mv HTAdapter73_mapped.sam STO18_mapped.sam
mv HTAdapter74_mapped.sam STO19_mapped.sam
mv HTAdapter75_mapped.sam STO20_mapped.sam
mv HTAdapter76_mapped.sam STO21_mapped.sam
mv HTAdapter77_mapped.sam SUM01_mapped.sam
mv HTAdapter78_mapped.sam SUM02_mapped.sam
mv HTAdapter79_mapped.sam SUM03_mapped.sam
mv HTAdapter80_mapped.sam SUM04_mapped.sam
mv HTAdapter81_mapped.sam SUM05_mapped.sam
mv HTAdapter82_mapped.sam SUM06_mapped.sam
mv HTAdapter83_mapped.sam SUM07_mapped.sam
mv HTAdapter84_mapped.sam SUM08_mapped.sam
mv HTAdapter85_mapped.sam SUM09_mapped.sam
mv HTAdapter86_mapped.sam SUM10_mapped.sam
mv HTAdapter87_mapped.sam SUM11_mapped.sam
mv HTAdapter88_mapped.sam SUM12_mapped.sam
mv HTAdapter89_mapped.sam SUM13_mapped.sam
mv HTAdapter90_mapped.sam SUM15_mapped.sam
mv HTAdapter91_mapped.sam SUM16_mapped.sam
mv HTAdapter92_mapped.sam SUM17_mapped.sam
mv HTAdapter93_mapped.sam SUM18_mapped.sam
mv HTAdapter94_mapped.sam SUM19_mapped.sam
mv HTAdapter95_mapped.sam SUM20_mapped.sam
mv HTAdapter96_mapped.sam SUM21_mapped.sam

## -- ADD READGROUPS -- ##
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_mapped.sam" O="ELF01_IDed.sam" RGPU="ELF" RGSM="ELF_54-38" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_mapped.sam" O="ELF02_IDed.sam" RGPU="ELF" RGSM="ELF_RA607" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_mapped.sam" O="ELF03_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_mapped.sam" O="ELF04_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_mapped.sam" O="ELF05_IDed.sam" RGPU="ELF" RGSM="ELF_558-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_mapped.sam" O="ELF06_IDed.sam" RGPU="ELF" RGSM="ELF_565-13" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_mapped.sam" O="ELF07_IDed.sam" RGPU="ELF" RGSM="ELF_5651-15" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_mapped.sam" O="ELF08_IDed.sam" RGPU="ELF" RGSM="ELF_568-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_mapped.sam" O="ELF09_IDed.sam" RGPU="ELF" RGSM="ELF_546-03" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_mapped.sam" O="ELF10_IDed.sam" RGPU="ELF" RGSM="ELF_547-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_mapped.sam" O="ELF11_IDed.sam" RGPU="ELF" RGSM="ELF_6050-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_mapped.sam" O="ELF12_IDed.sam" RGPU="ELF" RGSM="ELF_622-5" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_mapped.sam" O="ELF13_IDed.sam" RGPU="ELF" RGSM="ELF_636-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_mapped.sam" O="ELF14_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_mapped.sam" O="ELF15_IDed.sam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_mapped.sam" O="ELF16_IDed.sam" RGPU="ELF" RGSM="ELF_RP54" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_mapped.sam" O="MAH01_IDed.sam" RGPU="MAH" RGSM="MAH_35-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_mapped.sam" O="MAH02_IDed.sam" RGPU="MAH" RGSM="MAH_520-02" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_mapped.sam" O="MAH03_IDed.sam" RGPU="MAH" RGSM="MAH_539-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_mapped.sam" O="MAH04_IDed.sam" RGPU="MAH" RGSM="MAH_540-07" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_mapped.sam" O="MAH05_IDed.sam" RGPU="MAH" RGSM="MAH_541-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_mapped.sam" O="MAH06_IDed.sam" RGPU="MAH" RGSM="MAH_542-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_mapped.sam" O="MAH07_IDed.sam" RGPU="MAH" RGSM="MAH_543-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_mapped.sam" O="MAH08_IDed.sam" RGPU="MAH" RGSM="MAH_RA620" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_mapped.sam" O="MAH09_IDed.sam" RGPU="MAH" RGSM="MAH_RA641" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_mapped.sam" O="MAH10_IDed.sam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_mapped.sam" O="MAH11_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_mapped.sam" O="MAH12_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_mapped.sam" O="MAH13_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_mapped.sam" O="MAH14_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_mapped.sam" O="MAH15_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_mapped.sam" O="MAR01_IDed.sam" RGPU="MAR" RGSM="MAR_RP6771" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_mapped.sam" O="MAR02_IDed.sam" RGPU="MAR" RGSM="MAR_RP686" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_mapped.sam" O="MAR03_IDed.sam" RGPU="MAR" RGSM="MAR_368- 2" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_mapped.sam" O="MAR04_IDed.sam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_mapped.sam" O="MAR05_IDed.sam" RGPU="MAR" RGSM="MAR_RA379" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_mapped.sam" O="MAR06_IDed.sam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_mapped.sam" O="MAR07_IDed.sam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_mapped.sam" O="MAR08_IDed.sam" RGPU="MAR" RGSM="MAR_RA605" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_mapped.sam" O="MAR09_IDed.sam" RGPU="MAR" RGSM="MAR_RA630" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_mapped.sam" O="MAR10_IDed.sam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_mapped.sam" O="PAP01_IDed.sam" RGPU="PAP" RGSM="PAP_561-4" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_mapped.sam" O="PAP02_IDed.sam" RGPU="PAP" RGSM="PAP_562-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_mapped.sam" O="PAP03_IDed.sam" RGPU="PAP" RGSM="PAP_564-2" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_mapped.sam" O="PAP04_IDed.sam" RGPU="PAP" RGSM="PAP_569-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_mapped.sam" O="PAP05_IDed.sam" RGPU="PAP" RGSM="PAP_5010-05" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_mapped.sam" O="PAP06_IDed.sam" RGPU="PAP" RGSM="PAP_5021-02" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_mapped.sam" O="PAP07_IDed.sam" RGPU="PAP" RGSM="PAP_505-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_mapped.sam" O="PAP08_IDed.sam" RGPU="PAP" RGSM="PAP_516-02" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_mapped.sam" O="PAP09_IDed.sam" RGPU="PAP" RGSM="PAP_5241-01" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_mapped.sam" O="PAP10_IDed.sam" RGPU="PAP" RGSM="PAP_28-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_mapped.sam" O="PAP11_IDed.sam" RGPU="PAP" RGSM="PAP_RA434" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_mapped.sam" O="PAP12_IDed.sam" RGPU="PAP" RGSM="PAP_RA647" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_mapped.sam" O="PAP13_IDed.sam" RGPU="PAP" RGSM="PAP_RA648" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_mapped.sam" O="PAP14_IDed.sam" RGPU="PAP" RGSM="PAP_RA649" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_mapped.sam" O="PAP15_IDed.sam" RGPU="PAP" RGSM="PAP_RP16" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_mapped.sam" O="STO01_IDed.sam" RGPU="STO" RGSM="STO_RP697" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_mapped.sam" O="STO02_IDed.sam" RGPU="STO" RGSM="STO_RP698" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_mapped.sam" O="STO03_IDed.sam" RGPU="STO" RGSM="STO_RA530" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_mapped.sam" O="STO04_IDed.sam" RGPU="STO" RGSM="STO_RA549" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_mapped.sam" O="STO05_IDed.sam" RGPU="STO" RGSM="STO_RP22" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_mapped.sam" O="STO06_IDed.sam" RGPU="STO" RGSM="STO_271" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_mapped.sam" O="STO07_IDed.sam" RGPU="STO" RGSM="STO_275" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_mapped.sam" O="STO08_IDed.sam" RGPU="STO" RGSM="STO_43-113" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_mapped.sam" O="STO09_IDed.sam" RGPU="STO" RGSM="STO_43-115" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_mapped.sam" O="STO10_IDed.sam" RGPU="STO" RGSM="STO_43-125" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_mapped.sam" O="STO11_IDed.sam" RGPU="STO" RGSM="STO_43-76" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_mapped.sam" O="STO12_IDed.sam" RGPU="STO" RGSM="STO_43-77" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_mapped.sam" O="STO13_IDed.sam" RGPU="STO" RGSM="STO_43-78" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_mapped.sam" O="STO14_IDed.sam" RGPU="STO" RGSM="STO_43-74" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_mapped.sam" O="STO15_IDed.sam" RGPU="STO" RGSM="STO_43-75" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_mapped.sam" O="STO17_IDed.sam" RGPU="STO" RGSM="STO_33-10" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_mapped.sam" O="STO18_IDed.sam" RGPU="STO" RGSM="STO_560-2" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_mapped.sam" O="STO19_IDed.sam" RGPU="STO" RGSM="STO_528-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_mapped.sam" O="STO20_IDed.sam" RGPU="STO" RGSM="STO_559-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_mapped.sam" O="STO21_IDed.sam" RGPU="STO" RGSM="STO_33-272" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_mapped.sam" O="SUM01_IDed.sam" RGPU="SUM" RGSM="SUM_41-103" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_mapped.sam" O="SUM02_IDed.sam" RGPU="SUM" RGSM="SUM_41-104" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_mapped.sam" O="SUM03_IDed.sam" RGPU="SUM" RGSM="SUM_RA645" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_mapped.sam" O="SUM04_IDed.sam" RGPU="SUM" RGSM="SUM_RP31" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_mapped.sam" O="SUM05_IDed.sam" RGPU="SUM" RGSM="SUM_53-255" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_mapped.sam" O="SUM06_IDed.sam" RGPU="SUM" RGSM="SUM_53-252" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_mapped.sam" O="SUM07_IDed.sam" RGPU="SUM" RGSM="SUM_53-253" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_mapped.sam" O="SUM08_IDed.sam" RGPU="SUM" RGSM="SUM_53-256" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_mapped.sam" O="SUM09_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_mapped.sam" O="SUM10_IDed.sam" RGPU="SUM" RGSM="SUM_644-2" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_mapped.sam" O="SUM11_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_mapped.sam" O="SUM12_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_mapped.sam" O="SUM13_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_mapped.sam" O="SUM15_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_mapped.sam" O="SUM16_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_mapped.sam" O="SUM17_IDed.sam" RGPU="SUM" RGSM="SUM_624-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_mapped.sam" O="SUM18_IDed.sam" RGPU="SUM" RGSM="SUM_625-1" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_mapped.sam" O="SUM19_IDed.sam" RGPU="SUM" RGSM="SUM_626-3" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_mapped.sam" O="SUM20_IDed.sam" RGPU="SUM" RGSM="SUM_643-2" RGPL="illumina" RGLB="SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_mapped.sam" O="SUM21_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL="illumina" RGLB="SeqCap2012"

# RNA Seq
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629651_mapped.sam" O="SRR629651_IDed.sam" RGPU="SUM" RGSM="SAMN01823448" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629599_mapped.sam" O="SRR629599_IDed.sam" RGPU="SUM" RGSM="SAMN01823447" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629652_mapped.sam" O="SRR629652_IDed.sam" RGPU="SUM" RGSM="SAMN01823449" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629653_mapped.sam" O="SRR629653_IDed.sam" RGPU="SUM" RGSM="SAMN01823450" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629654_mapped.sam" O="SRR629654_IDed.sam" RGPU="SUM" RGSM="SAMN01823451" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629655_mapped.sam" O="SRR629655_IDed.sam" RGPU="SUM" RGSM="SAMN01823452" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629656_mapped.sam" O="SRR629656_IDed.sam" RGPU="SUM" RGSM="SAMN01823453" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629657_mapped.sam" O="SRR629657_IDed.sam" RGPU="SUM" RGSM="SAMN01823454" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629658_mapped.sam" O="SRR629658_IDed.sam" RGPU="SUM" RGSM="SAMN01823455" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629659_mapped.sam" O="SRR629659_IDed.sam" RGPU="SUM" RGSM="SAMN01823456" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629661_mapped.sam" O="SRR629661_IDed.sam" RGPU="SUM" RGSM="SAMN01823457" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629660_mapped.sam" O="SRR629660_IDed.sam" RGPU="SUM" RGSM="SAMN01823458" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629662_mapped.sam" O="SRR629662_IDed.sam" RGPU="SUM" RGSM="SAMN01823459" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629664_mapped.sam" O="SRR629664_IDed.sam" RGPU="SUM" RGSM="SAMN01823460" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629663_mapped.sam" O="SRR629663_IDed.sam" RGPU="SUM" RGSM="SAMN01823461" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629665_mapped.sam" O="SRR629665_IDed.sam" RGPU="SUM" RGSM="SAMN01823462" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629666_mapped.sam" O="SRR629666_IDed.sam" RGPU="SUM" RGSM="SAMN01823463" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629667_mapped.sam" O="SRR629667_IDed.sam" RGPU="SUM" RGSM="SAMN01823464" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629668_mapped.sam" O="SRR629668_IDed.sam" RGPU="SUM" RGSM="SAMN01823465" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629669_mapped.sam" O="SRR629669_IDed.sam" RGPU="SUM" RGSM="SAMN01823466" RGPL="illumina" RGLB="RNASeq2012" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497737_mapped.sam" O="SRR497737_IDed.sam" RGPU="SUM" RGSM="SAMN00996377" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497739_mapped.sam" O="SRR497739_IDed.sam" RGPU="SUM" RGSM="SAMN00996379" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497740_mapped.sam" O="SRR497740_IDed.sam" RGPU="SUM" RGSM="SAMN00996380" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497741_mapped.sam" O="SRR497741_IDed.sam" RGPU="SUM" RGSM="SAMN00996381" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497742_mapped.sam" O="SRR497742_IDed.sam" RGPU="SUM" RGSM="SAMN00996382" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497743_mapped.sam" O="SRR497743_IDed.sam" RGPU="SUM" RGSM="SAMN00996383" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497744_mapped.sam" O="SRR497744_IDed.sam" RGPU="SUM" RGSM="SAMN00996384" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497745_mapped.sam" O="SRR497745_IDed.sam" RGPU="SUM" RGSM="SAMN00996385" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497746_mapped.sam" O="SRR497746_IDed.sam" RGPU="SUM" RGSM="SAMN00996386" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497747_mapped.sam" O="SRR497747_IDed.sam" RGPU="SUM" RGSM="SAMN00996387" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497748_mapped.sam" O="SRR497748_IDed.sam" RGPU="SUM" RGSM="SAMN00996388" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497749_mapped.sam" O="SRR497749_IDed.sam" RGPU="SUM" RGSM="SAMN00996389" RGPL="illumina" RGLB="RNASeq2008" 
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497738_mapped.sam" O="SRR497738_IDed.sam" RGPU="SUM" RGSM="SAMN00996378" RGPL="illumina" RGLB="RNASeq2008" 

# -- PREPARE FOR SNP CALLING -- ##
# HybPiper
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA
# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# index ref
samtools faidx /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fna
java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
    R= /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fna \
    O= /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.dict 
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA
# Create gvcf map
touch /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/cohort.sample_map
while read i;
do
# convert .sam to .bam & sort ###
samtools view -@ 2 -bS /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/"$i"_sorted.bam
### remove duplicates ###
java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
    INPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/"$i"_sorted.bam \
    OUTPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam \
    METRICS_FILE=DuplicationMetrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORTED=TRUE \
    REMOVE_DUPLICATES=TRUE
#index sorted bam
samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam
### ---------------- Variant Calling ------------------- ###
/tools/gatk-4.0.10.1/gatk --java-options "-Xmx16g" HaplotypeCaller \
	-R /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fna
	-I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam \
	-O /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_rawVariants.g.vcf
	-ERC GVCF
echo "$i\t$i_rawVariants.g.vcf" >> /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/cohort.sample_map 
	
## Calculate Mapping Stats ##
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/ 
# tally mapped reads & calcuate the stats
samtools idxstats /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/StatsDNA/"$i"_counts.txt
samtools flagstat /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/StatsDNA/"$i"_stats.txt
samtools depth /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/StatsDNA/"$i"_depth.txt
done<samList
 
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/StatsDNA
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList
 
# Generate file with all depth information
cat *_depth.txt > DNA_depth.txt
for f in DNA_depth.txt
do
sed -i "s/$/\tDNA/" $f; done

#  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
 #  # make results directory & move results
 #  mkdir -p /home/rlk0015/SeqCap/May2020/stats
 #  mkdir -p /home/rlk0015/SeqCap/May2020/counts
 #  mkdir -p /home/rlk0015/SeqCap/May2020/avgDepth
 #  
 #  # cp results to respective directories
 #  cp *stats.txt /home/rlk0015/SeqCap/May2020/stats
 #  cp *counts.txt /home/rlk0015/SeqCap/May2020/counts
 #  cp *depth.txt /home/rlk0015/SeqCap/May2020/avgDepth
 #  *****************************************************************************************************

### --------------------- Continue Variant Calling ----------------------- ###

# Combine GVCFs
/tools/gatk-4.0.10.1/gatk --java-options "-Xmx16g" GenomicsDBImport \
	--sample_name_map /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/cohort.sample_map
	--genomicsdb-workspace-path SNP_database
	--reader-threads 5
# HERE 
# ### move to GATK directory ###
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/
# ### merge .bam files ###
# samtools merge -f /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/dupsRemoved.bam /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/*_dupsRemoved.bam
# ### index the merged .bam ###
#  samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/dupsRemoved.bam
#  # realign indels
#  java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T IndelRealigner \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/dupsRemoved.bam \
#      -targetIntervals /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsCalled.intervals \
#      -LOD 3.0 \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsRealigned.bam
#  ## -- CALL SNPS -- ##
#  java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T UnifiedGenotyper \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsRealigned.bam \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/calledSNPs.vcf \
#      -gt_mode DISCOVERY \
#      -ploidy 2 \
#      -stand_call_conf 30 \
#      -stand_emit_conf 10 \
#      -rf BadCigar
#  # annotate variants
#  java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T VariantAnnotator \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsRealigned.bam \
#      -G StandardAnnotation \
#      -V:variant,VCF /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/calledSNPs.vcf \
#      -XA SnpEff \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/annotatedSNPs.vcf \
#      -rf BadCigar
#  # annotate indels
#  java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T UnifiedGenotyper \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsRealigned.bam \
#      -gt_mode DISCOVERY \
#      -glm INDEL \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/annotatedIndels.vcf \
#      -stand_call_conf 30 \
#      -stand_emit_conf 10 \
#      -rf BadCigar
#  # mask indels
#  java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T VariantFiltration \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/calledSNPs.vcf \
#      --mask /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/annotatedIndels.vcf \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/SNPsMaskedIndels.vcf/ \
#      -rf BadCigar
#  # restrict to high-quality variant calls
#  cat /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/qualitySNPs.vcf
#  # read-backed phasing
#  java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T ReadBackedPhasing \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/indelsRealigned.bam \
#      --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/qualitySNPs.vcf \
#      -L /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/qualitySNPs.vcf \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/phasedSNPs.vcf \
#      --phaseQualityThresh 20.0 \
#      -rf BadCigar
#  
#  ## Make Sample List from VCF ##
#  cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA
#  bcftools query -l phasedSNPs.vcf > VcfSampleList
#  
#  while read i;
#  do
#  # VCF for each sample
#  java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T SelectVariants \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/phasedSNPs.vcf \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_phasedSNPs.vcf \
#      -sn "$i" \
#      -rf BadCigar
#  # make SNPs table
#  java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#      -T VariantsToTable \
#      -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/"$i"_phasedSNPs.vcf \
#      -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
#      -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/HybPiperSNPTables/"$i"_tableSNPs.txt \
#      -rf BadCigar
#  # Add phased SNPs to reference and filter
#  python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
#      /scratch/rlk0015/Telag/May2020/WorkingDirectory/GATKDNA/HybPiperContigs.fa \
#      /scratch/rlk0015/Telag/May2020/WorkingDirectory/HybPiperSNPTables/"$i"_tableSNPs.txt \
#      /scratch/rlk0015/Telag/May2020/WorkingDirectory/HybPiperSequenceTables/"$i"_tableSequences.txt \
#      1
#  done<VcfSampleList
# 
# # ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# # ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# # ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# # ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# # ºººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº #
# 
# 
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ TRANSCRIPTS BLOCK ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Following is complete &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/Transcripts
# # && # make .sam files list for Samtools processing
# # && ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# # && # make result directories
# mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK
# mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/Sequences
# # && mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsSNPTables
# # && mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsSequenceTables
# # && # copy reference for SNP calling
# cp /home/rlk0015/SeqCap/code/References/Transcripts.fa /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa
# # && # index ref
# # && samtools faidx /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa
# # && java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
# # &&     R= /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     O= /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.dict
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/Transcripts
# # && while read i;
# # && do
# # && # convert .sam to .bam & sort ###
# # && samtools view -@ 2 -bS /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"_sorted.bam
# # && ### remove duplicates ###
# # && java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
# # &&     INPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsDNA/Transcripts/"$i"_sorted.bam \
# # &&     OUTPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam \
# # &&     METRICS_FILE=DuplicationMetrics \
# # &&     CREATE_INDEX=true \
# # &&     VALIDATION_STRINGENCY=SILENT \
# # &&     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
# # &&     ASSUME_SORTED=TRUE \
# # &&     REMOVE_DUPLICATES=TRUE
# # && #index sorted bam
# # && samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam
# # && 
# # && ## Calculate Stats ##
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory_ReducedTranscriptome/mappedReadsDNA/Transcripts
# # && # make stats folder
# # && mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsStats
# # && # tally mapped reads & calcuate the stats
# # && samtools idxstats /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsStats/"$i"_counts.txt
# # && samtools flagstat /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsStats/"$i"_stats.txt
# # && samtools depth /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsStats/"$i"_depth.txt
# # && done<samList
# # && 
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsStats
# # && ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# # && # Add individual name to each line in depth file
# # && while read i
# # && do
# # && for f in "$i"_depth.txt
# # && do
# # && sed -i "s/$/\t$i/" $f; done
# # && done<depthList
# # && 
# # && # Generate file with all depth information
# # && cat *_depth.txt > Transcripts_depth.txt
# # && for f in Transcripts_depth.txt
# # && do
# # && sed -i "s/$/\tTranscripts/" $f; done
# # && 
# # && # Create file with avg depth per exon using Randy's python script
# # && #!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
# # && python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py Transcripts_depth.txt Transcripts_avgDepth.txt
# # && 
# # && #  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
# # && #  # make results directory & move results
# # && #  mkdir -p /home/rlk0015/SeqCap/May2020/stats
# # && #  mkdir -p /home/rlk0015/SeqCap/May2020/counts
# # && #  mkdir -p /home/rlk0015/SeqCap/May2020/avgDepth
# # && #  
# # && #  # cp results to respective directories
# # && #  cp *stats.txt /home/rlk0015/SeqCap/May2020/stats
# # && #  cp *counts.txt /home/rlk0015/SeqCap/May2020/counts
# # && #  cp *depth.txt /home/rlk0015/SeqCap/May2020/avgDepth
# # && #  *****************************************************************************************************
# # && 
# # && 
# # && ### move to GATK directory ###
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/
# # && ### merge .bam files ###
# # && samtools merge -f /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/*_dupsRemoved.bam
# # && ### index the merged .bam ###
# # && samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/dupsRemoved.bam
# # && ## -- CALL SNPS -- ##
# /tools/gatk-4.1.2.0/gatk --java-options "-Xmx16g" HaplotypeCaller\
#   -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
#   -I dupsRemoved.bam \
#   -stand_emit_conf 10 \
#   -stand_call_conf 30 \
#   -O rawVariants.vcf 
#   # The following read filters are applied automatically: HCMappingQualityFilter, MalformedReadFilter, BadCigarFilter, UnmappedReadFilter, NotPrimaryAlignmentFilter, FailsVendorQualityCheckFilter, DuplicateReadFilter, MappingQualityUnavailableFilter
# # annotate variants
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T VariantAnnotator \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
# # &&     -G StandardAnnotation \
# # &&     -V:variant,VCF /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/calledSNPs.vcf \
# # &&     -XA SnpEff \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/annotatedSNPs.vcf \
# # &&     -rf BadCigar
# # && # annotate indels
# # && java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T UnifiedGenotyper \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
# # &&     -gt_mode DISCOVERY \
# # &&     -glm INDEL \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/annotatedIndels.vcf \
# # &&     -stand_call_conf 30 \
# # &&     -stand_emit_conf 10 \
# # &&     -rf BadCigar
# # && # mask indels
# # && java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T VariantFiltration \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/calledSNPs.vcf \
# # &&     --mask /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/annotatedIndels.vcf \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/SNPsMaskedIndels.vcf/ \
# # &&     -rf BadCigar
# # && # restrict to high-quality variant calls
# # && cat /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf
# # && # read-backed phasing
# # && java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T ReadBackedPhasing \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/indelsRealigned.bam \
# # &&     --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf \
# # &&     -L /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/qualitySNPs.vcf \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/phasedSNPs.vcf \
# # &&     --phaseQualityThresh 20.0 \
# # &&     -rf BadCigar
# # && 
# # && ## Make Sample List from VCF ##
# # && cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK
# # && bcftools query -l phasedSNPs.vcf > VcfSampleList
# # && 
# # && while read i;
# # && do
# # && # VCF for each sample
# # && java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T SelectVariants \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/phasedSNPs.vcf \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_phasedSNPs.vcf \
# # &&     -sn "$i" \
# # &&     -rf BadCigar
# # && # make SNPs table
# # && java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
# # &&     -T VariantsToTable \
# # &&     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/"$i"_phasedSNPs.vcf \
# # &&     -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
# # &&     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsSNPTables/"$i"_tableSNPs.txt \
# # &&     -rf BadCigar
# # && # Add phased SNPs to reference and filter
# # && python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
# # &&     /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsGATK/Transcripts.fa \
# # &&     /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsSNPTables/"$i"_tableSNPs.txt \
# # &&     /scratch/rlk0015/Telag/May2020/WorkingDirectory/TranscriptsSequenceTables/"$i"_tableSequences.txt \
# # &&     1
# # && done<VcfSampleList
# # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& End && block &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# 
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ RNA SEQ BLOCK ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA
# # +++++++++++++++++++++++++ Already done +++++++++++++++++++++++++++++
# # ++ # make .sam files list for Samtools processing
# # ++ ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# # ++ # make result directories
# # ++ mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK
# # ++ mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNASNPTables
# # ++ mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNASequenceTables
# # ++ # copy reference for SNP calling
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cp /home/rlk0015/SeqCap/code/References/TranscriptomePlusTranscripts.fa /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa
# # index ref
# samtools faidx /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa
# java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
#     R= /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     O= /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.dict
# # =========================== Already done ===========================
# # == cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA
# # == while read i;
# # == do
# # == # convert .sam to .bam & sort ###
# # == samtools view -@ 2 -bS /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA/"$i"*_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam
# # == ### remove duplicates ###
# # == java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
# # ==     INPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam \
# # ==     OUTPUT=/scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam \
# # ==     METRICS_FILE=DuplicationMetrics \
# # ==     CREATE_INDEX=true \
# # ==     VALIDATION_STRINGENCY=SILENT \
# # ==     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
# # ==     ASSUME_SORTED=TRUE \
# # ==     REMOVE_DUPLICATES=TRUE
# # == #index sorted bam
# # == samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam
# # == 
# # == ## Calculate Stats ##
# # == cd /scratch/rlk0015/Telag/May2020/WorkingDirectory_ReducedTranscriptome/mappedReadsRNA
# # == # make stats folder
# # == mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAStats
# # == # tally mapped reads & calcuate the stats
# # == samtools idxstats /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAStats/"$i"_counts.txt
# # == samtools flagstat /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAStats/"$i"_stats.txt
# # == samtools depth /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAStats/"$i"_depth.txt
# # == done<samList
# # == 
# # == cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAStats
# # == ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# # == # Add individual name to each line in depth file
# # == while read i
# # == do
# # == for f in "$i"_depth.txt
# # == do
# # == sed -i "s/$/\t$i/" $f; done
# # == done<depthList
# # == 
# # == # Generate file with all depth information
# # == cat *_depth.txt > RNA_depth.txt
# # == for f in RNA_depth.txt
# # == do
# # == sed -i "s/$/\tRNA/" $f; done
# # == 
# # == # Create file with avg depth per exon using Randy's python script
# # == #!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
# # == python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py RNA_depth.txt RNA_avgDepth.txt
# # == 
# # == #  ************************ Commented out because I don't want to worry about this right now- needs to be done eventually *******
# # == #  # make results directory & move results
# # == #  mkdir -p /home/rlk0015/SeqCap/May2020/stats
# # == #  mkdir -p /home/rlk0015/SeqCap/May2020/counts
# # == #  mkdir -p /home/rlk0015/SeqCap/May2020/avgDepth
# # == #  
# # == #  # cp results to respective directories
# # == #  cp *stats.txt /home/rlk0015/SeqCap/May2020/stats
# # == #  cp *counts.txt /home/rlk0015/SeqCap/May2020/counts
# # == #  cp *depth.txt /home/rlk0015/SeqCap/May2020/avgDepth
# # == #  *****************************************************************************************************
# # == =================================================================================== 
# 
# 
# ### move to GATK directory ###
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/
# # .......................... Already done ..............................
# # .. ### merge .bam files ###
# # .. samtools merge -f /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/dupsRemoved.bam /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/*_dupsRemoved.bam
# # .. ### index the merged .bam ###
# # .. samtools index /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/dupsRemoved.bam
# # call indels
# # ......................................................................
# java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T RealignerTargetCreator \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/dupsRemoved.bam \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsCalled.intervals
# # realign indels
# java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T IndelRealigner \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/dupsRemoved.bam \
#     -targetIntervals /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsCalled.intervals \
#     -LOD 3.0 \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsRealigned.bam
# ## -- CALL SNPS -- ##
# java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T UnifiedGenotyper \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsRealigned.bam \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/calledSNPs.vcf \
#     -gt_mode DISCOVERY \
#     -ploidy 2 \
#     -stand_call_conf 30 \
#     -stand_emit_conf 10 \
#     -rf BadCigar
# # annotate variants
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T VariantAnnotator \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsRealigned.bam \
#     -G StandardAnnotation \
#     -V:variant,VCF /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/calledSNPs.vcf \
#     -XA SnpEff \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/annotatedSNPs.vcf \
#     -rf BadCigar
# # annotate indels
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T UnifiedGenotyper \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsRealigned.bam \
#     -gt_mode DISCOVERY \
#     -glm INDEL \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/annotatedIndels.vcf \
#     -stand_call_conf 30 \
#     -stand_emit_conf 10 \
#     -rf BadCigar
# # mask indels
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T VariantFiltration \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/calledSNPs.vcf \
#     --mask /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/annotatedIndels.vcf \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/SNPsMaskedIndels.vcf/ \
#     -rf BadCigar
# # restrict to high-quality variant calls
# cat /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/qualitySNPs.vcf
# # read-backed phasing
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T ReadBackedPhasing \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -I /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/indelsRealigned.bam \
#     --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/qualitySNPs.vcf \
#     -L /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/qualitySNPs.vcf \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/phasedSNPs.vcf \
#     --phaseQualityThresh 20.0 \
#     -rf BadCigar
# 
# ## Make Sample List from VCF ##
# cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK
# bcftools query -l phasedSNPs.vcf > VcfSampleList
# 
# while read i;
# do
# # VCF for each sample
# java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T SelectVariants \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     --variant /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/phasedSNPs.vcf \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_phasedSNPs.vcf \
#     -sn "$i" \
#     -rf BadCigar
# # make SNPs table
# java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
#     -T VariantsToTable \
#     -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     -V /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/"$i"_phasedSNPs.vcf \
#     -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
#     -o /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNASNPTables/"$i"_tableSNPs.txt \
#     -rf BadCigar
# # Add phased SNPs to reference and filter
# python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
#     /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNAGATK/Transcriptome.fa \
#     /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNASNPTables/"$i"_tableSNPs.txt \
#     /scratch/rlk0015/Telag/May2020/WorkingDirectory/RNASequenceTables/"$i"_tableSequences.txt \
#     1
# done<VcfSampleList
