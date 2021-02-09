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
#PBS -l nodes=3:ppn=10,walltime=05:00:00:00
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
module load bcftools/1.3.2
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

#+   ++++++++++++++++++++++++ Done June 29: 1994303 & 1994412 ++++++++++++++++++++++++   +#
#+   #create working directory ###
#+   mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory 
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
#+   # DNA
#+   mkdir -p $WorkingDirectory/rawReadsDNA
#+   mkdir -p $WorkingDirectory/cleanReadsDNA
#+   mkdir -p $WorkingDirectory/mappedReadsDNA
#+   mkdir -p $WorkingDirectory/StatsDNA
#+   mkdir -p $WorkingDirectory/GATKDNA
#+   mkdir -p $WorkingDirectory/SNPTablesDNA
#+   mkdir -p $WorkingDirectory/SequenceTablesDNA
#+   # RNA
#+   mkdir -p $WorkingDirectory/rawReadsRNA
#+   mkdir -p $WorkingDirectory/cleanReadsRNA
#+   mkdir -p $WorkingDirectory/mappedReadsRNA
#+   mkdir -p $WorkingDirectory/StatsRNA
#+   mkdir -p $WorkingDirectory/GATKRNA
#+   mkdir -p $WorkingDirectory/SNPTablesRNA
#+   mkdir -p $WorkingDirectory/SequenceTablesRNA
#+   # Utilities
#+   mkdir -p $WorkingDirectory/References
#+   echo "RLK_report: directory created: $WorkingDirectory with rawReads and cleanReads sub directories"
#+   
#+   #+   ++++++++++++++++++++++++ Done June 27: 1994197 ++++++++++++++++++++++++   +#
#+   #+   # ---------------------------
#+   #+   # Copy raw DNA reads over
#+   #+   # ---------------------------
#+   #+   cd /home/shared/tss0019_lab/SeqCap_GarterSnake2012/
#+   #+   for i in {1..96}
#+   #+   do
#+   #+     cp Sample_HTAdapter"$i"/*.fastq.gz $WorkingDirectory/rawReadsDNA
#+   #+   done
#+   #+   echo "RLK_report: SEQ-CAP RAW READ COPY COMPLETE"
#+   #+   
#+   #+   ### perform initial quality check ###
#+   #+   cd $WorkingDirectory/rawReadsDNA
#+   #+   ls *.fastq.gz | time parallel -j+0 --eta 'fastqc {}'
#+   #+   multiqc .
#+   #+   echo "RLK_report: SEQ-CAP RAW READ QUALITY CHECK COMPLETE"
#+   #+   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   +#
#+   
#+   ### copy over adapter file ###
#+   cp /home/rlk0015/SeqCap/code/References/adapters.fa $WorkingDirectory/rawReadsDNA
#+   echo "RLK_report: ADAPTERS COPY COMPLETE"
#+   
#+   cd $WorkingDirectory/rawReadsDNA
#+   ### paired-end trimming ###
#+   ls | grep "fastq.gz" | cut -d "_" -f 1,2 | sort | uniq > PE_TrimmList
#+   while read i
#+   do
#+   java -jar /tools/trimmomatic-0.36/trimmomatic-0.36.jar \
#+     PE \
#+     -threads 6 \
#+     -phred33 \
#+     $WorkingDirectory/rawReadsDNA/"$i"*_R1*.fastq.gz \
#+     $WorkingDirectory/rawReadsDNA/"$i"*_R2*.fastq.gz \
#+     $WorkingDirectory/cleanReadsDNA/"$i"_R1_paired.fastq.gz \ 
#+     $WorkingDirectory/cleanReadsDNA/"$i"_R1_unpaired.fastq.gz \
#+     $WorkingDirectory/cleanReadsDNA/"$i"_R2_paired.fastq.gz \
#+     $WorkingDirectory/cleanReadsDNA/"$i"_R2_unpaired.fastq.gz \
#+     ILLUMINACLIP:adapters.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:30 MINLEN:36
#+   done<PE_TrimmList
#+   echo "RLK_report: SEQ-CAP READ TRIMMING COMPLETE"
#+   ### perform second quality check ###
#+   cd $WorkingDirectory/cleanReadsDNA/
#+   ls *paired.fastq.gz | grep "paired.fastq.gz" | time parallel -j+0 --eta 'fastqc {}'
#+   multiqc .
#+   echo "RLK_report: SEQ-CAP CLEAN READ QUALITY CHECK COMPLETE"
#+   
#+   # ---------------------------
#+   # Copy clean RNA reads over
#+   # ---------------------------
#+   cd /scratch/GarterSnake/RNAseq_2012Data/CleanedData
#+   ls SRR*.fastq | parallel -j+0 --eta 'cp {} '$WorkingDirectory'/cleanReadsRNA'
#+   cd /scratch/GarterSnake/RNAseq_2008Data/CleanedData
#+   ls SRR*.fastq | parallel -j+0 --eta 'cp {} '$WorkingDirectory'/cleanReadsRNA'
#+   echo "RLK_report: RNA-SEQ CLEAN READ COPY COMPLETE"
#+   
#+   # -------------- Already performed- see details below --------------
#+   # # Tonia's criteria for trimming RNA-Seq reads:
#+   # ########## 2008 Data Trimmomatic #############
#+   # java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar SE -phred33 "$i"_1.fastq  /scratch/GarterSnake/RNAseq_2008Data/CleanedData/"$i"_cleaned.fastq  LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:36 2012 samples are PE 100 bp reads.
#+   # ########## 2012 Trimmomatic #############
#+   # java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar PE  -phred33 "$i"_1.fastq "$i"_2.fastq "$i"_1_paired.fastq "$i"_1_unpaired.fastq "$i"_2_paired.fastq "$i"_2_unpaired.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
#+   # ------------------------------------------------------------------
#+   
#+   # -------------------------------------
#+   # Quality check on RNA-Seq cleaned data 
#+   # -------------------------------------
#+   cd $WorkingDirectory/cleanReadsRNA
#+   ls *.fastq | grep "fastq" | time parallel -j+0 --eta 'fastqc {}'
#+   multiqc . 
#+   echo "RLK_report: RNA-SEQ CLEAN READ QUALITY CHECK COMPLETE"
#+   
#+   
#+   # ------------------------
#+   # Map Seq-Cap data to ref
#+   # ------------------------
#+   ### ••••••••••••••••••••••• This section is complete ••••••••••••••••••••• ###
#+   # ### copy the reference (T. elegans genome) to References directory ###
#+   # cp /home/rlk0015/SeqCap/code/References/T_elegans_genome/latest_assembly_versions/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.fna.gz $WorkingDirectory/References/TelagGenome.fasta.gz
#+   # cd $WorkingDirectory/References
#+   # gunzip TelagGenome.fasta.gz
#+   # echo "RLK_report: REFERENCE GENOME COPY COMPLETE"
#+   # index T. elegans genome
#+   # cd $WorkingDirectory/References
#+   # bwa index -p TelagGenome -a bwtsw TelagGenome.fasta
#+   # echo "RLK_report: REFERENCE GENOME BWA INDEX COMPLETE"
#+   ### •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• ###
#+   
#+   # create list with each paired individual
#+   cd $WorkingDirectory/cleanReadsDNA
#+   ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > \
#+   	$WorkingDirectory/mappedReadsDNA/pairedMapList
#+   echo "pairedMapList created"
#+   cd $WorkingDirectory/mappedReadsDNA
#+   # while loop through the names in pairedMapList
#+   while read i
#+   do
#+   ### map to T. elegans genome ###
#+   bwa mem -t 4 -M $WorkingDirectory/References/TelagGenome \
#+   	$WorkingDirectory/cleanReadsDNA/"$i"_*R1_paired.fastq.gz \
#+   	$WorkingDirectory/cleanReadsDNA/"$i"_*R2_paired.fastq.gz > \
#+   	$WorkingDirectory/mappedReadsDNA/"$i"_mapped.sam
#+   done<$WorkingDirectory/mappedReadsDNA/pairedMapList
#+   echo "mapped DNA reads to genome"
#+   
#+   ### ***  -------------------  MAP RNA TO REFERENCE  ------------------ *** ###
#+   # create list with each paired individual
#+   cd $WorkingDirectory/cleanReadsRNA
#+   ls | grep "_cleaned.fastq" | cut -d "_" -f 1 | sort | uniq > $WorkingDirectory/mappedReadsRNA/singleEndMapList
#+   ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > $WorkingDirectory/mappedReadsRNA/pairedEndMapList
#+   cd $WorkingDirectory/mappedReadsRNA
#+   # while loop through the names in pairedEndMapList
#+   while read i
#+   do
#+   ### PAIRED-END MAPPING ###
#+   bwa mem -t 4 -M $WorkingDirectory/References/TelagGenome \
#+     $WorkingDirectory/cleanReadsRNA/"$i"*_1_paired.fastq \
#+     $WorkingDirectory/cleanReadsRNA/"$i"*_2_paired.fastq > \
#+     $WorkingDirectory/mappedReadsRNA/"$i"_mapped.sam
#+   done<$WorkingDirectory/mappedReadsRNA/pairedEndMapList
#+   echo "RLK_report: paired map complete"
#+   ### SINGLE-END MAPPING ###
#+   while read i
#+   do
#+   bwa mem -t 4 -M $WorkingDirectory/References/TelagGenome \
#+     $WorkingDirectory/cleanReadsRNA/"$i"*_cleaned.fastq > \
#+     $WorkingDirectory/mappedReadsRNA/"$i"_mapped.sam
#+   done<$WorkingDirectory/mappedReadsRNA/singleEndMapList
#+   
#+   ### change names to match population sampling ###
#+   cd $WorkingDirectory/mappedReadsDNA/
#+   mv HTAdapter1_mapped.sam ELF01_mapped.sam
#+   mv HTAdapter2_mapped.sam ELF02_mapped.sam
#+   mv HTAdapter3_mapped.sam ELF03_mapped.sam
#+   mv HTAdapter4_mapped.sam ELF04_mapped.sam
#+   mv HTAdapter5_mapped.sam ELF05_mapped.sam
#+   mv HTAdapter6_mapped.sam ELF06_mapped.sam
#+   mv HTAdapter7_mapped.sam ELF07_mapped.sam
#+   mv HTAdapter8_mapped.sam ELF08_mapped.sam
#+   mv HTAdapter9_mapped.sam ELF09_mapped.sam
#+   mv HTAdapter10_mapped.sam ELF10_mapped.sam
#+   mv HTAdapter11_mapped.sam ELF11_mapped.sam
#+   mv HTAdapter12_mapped.sam ELF12_mapped.sam
#+   mv HTAdapter13_mapped.sam ELF13_mapped.sam
#+   mv HTAdapter14_mapped.sam ELF14_mapped.sam
#+   mv HTAdapter15_mapped.sam ELF15_mapped.sam
#+   mv HTAdapter16_mapped.sam ELF16_mapped.sam
#+   mv HTAdapter17_mapped.sam MAH01_mapped.sam
#+   mv HTAdapter18_mapped.sam MAH02_mapped.sam
#+   mv HTAdapter19_mapped.sam MAH03_mapped.sam
#+   mv HTAdapter20_mapped.sam MAH04_mapped.sam
#+   mv HTAdapter21_mapped.sam MAH05_mapped.sam
#+   mv HTAdapter22_mapped.sam MAH06_mapped.sam
#+   mv HTAdapter23_mapped.sam MAH07_mapped.sam
#+   mv HTAdapter24_mapped.sam MAH08_mapped.sam
#+   mv HTAdapter25_mapped.sam MAH09_mapped.sam
#+   mv HTAdapter26_mapped.sam MAH10_mapped.sam
#+   mv HTAdapter27_mapped.sam MAH11_mapped.sam
#+   mv HTAdapter28_mapped.sam MAH12_mapped.sam
#+   mv HTAdapter29_mapped.sam MAH13_mapped.sam
#+   mv HTAdapter30_mapped.sam MAH14_mapped.sam
#+   mv HTAdapter31_mapped.sam MAH15_mapped.sam
#+   mv HTAdapter32_mapped.sam MAR01_mapped.sam
#+   mv HTAdapter33_mapped.sam MAR02_mapped.sam
#+   mv HTAdapter34_mapped.sam MAR03_mapped.sam
#+   mv HTAdapter35_mapped.sam MAR04_mapped.sam
#+   mv HTAdapter36_mapped.sam MAR05_mapped.sam
#+   mv HTAdapter37_mapped.sam MAR06_mapped.sam
#+   mv HTAdapter38_mapped.sam MAR07_mapped.sam
#+   mv HTAdapter39_mapped.sam MAR08_mapped.sam
#+   mv HTAdapter40_mapped.sam MAR09_mapped.sam
#+   mv HTAdapter41_mapped.sam MAR10_mapped.sam
#+   mv HTAdapter42_mapped.sam PAP01_mapped.sam
#+   mv HTAdapter43_mapped.sam PAP02_mapped.sam
#+   mv HTAdapter44_mapped.sam PAP03_mapped.sam
#+   mv HTAdapter45_mapped.sam PAP04_mapped.sam
#+   mv HTAdapter46_mapped.sam PAP05_mapped.sam
#+   mv HTAdapter47_mapped.sam PAP06_mapped.sam
#+   mv HTAdapter48_mapped.sam PAP07_mapped.sam
#+   mv HTAdapter49_mapped.sam PAP08_mapped.sam
#+   mv HTAdapter50_mapped.sam PAP09_mapped.sam
#+   mv HTAdapter51_mapped.sam PAP10_mapped.sam
#+   mv HTAdapter52_mapped.sam PAP11_mapped.sam
#+   mv HTAdapter53_mapped.sam PAP12_mapped.sam
#+   mv HTAdapter54_mapped.sam PAP13_mapped.sam
#+   mv HTAdapter55_mapped.sam PAP14_mapped.sam
#+   mv HTAdapter56_mapped.sam PAP15_mapped.sam
#+   mv HTAdapter57_mapped.sam STO01_mapped.sam
#+   mv HTAdapter58_mapped.sam STO02_mapped.sam
#+   mv HTAdapter59_mapped.sam STO03_mapped.sam
#+   mv HTAdapter60_mapped.sam STO04_mapped.sam
#+   mv HTAdapter61_mapped.sam STO05_mapped.sam
#+   mv HTAdapter62_mapped.sam STO06_mapped.sam
#+   mv HTAdapter63_mapped.sam STO07_mapped.sam
#+   mv HTAdapter64_mapped.sam STO08_mapped.sam
#+   mv HTAdapter65_mapped.sam STO09_mapped.sam
#+   mv HTAdapter66_mapped.sam STO10_mapped.sam
#+   mv HTAdapter67_mapped.sam STO11_mapped.sam
#+   mv HTAdapter68_mapped.sam STO12_mapped.sam
#+   mv HTAdapter69_mapped.sam STO13_mapped.sam
#+   mv HTAdapter70_mapped.sam STO14_mapped.sam
#+   mv HTAdapter71_mapped.sam STO15_mapped.sam
#+   mv HTAdapter72_mapped.sam STO17_mapped.sam
#+   mv HTAdapter73_mapped.sam STO18_mapped.sam
#+   mv HTAdapter74_mapped.sam STO19_mapped.sam
#+   mv HTAdapter75_mapped.sam STO20_mapped.sam
#+   mv HTAdapter76_mapped.sam STO21_mapped.sam
#+   mv HTAdapter77_mapped.sam SUM01_mapped.sam
#+   mv HTAdapter78_mapped.sam SUM02_mapped.sam
#+   mv HTAdapter79_mapped.sam SUM03_mapped.sam
#+   mv HTAdapter80_mapped.sam SUM04_mapped.sam
#+   mv HTAdapter81_mapped.sam SUM05_mapped.sam
#+   mv HTAdapter82_mapped.sam SUM06_mapped.sam
#+   mv HTAdapter83_mapped.sam SUM07_mapped.sam
#+   mv HTAdapter84_mapped.sam SUM08_mapped.sam
#+   mv HTAdapter85_mapped.sam SUM09_mapped.sam
#+   mv HTAdapter86_mapped.sam SUM10_mapped.sam
#+   mv HTAdapter87_mapped.sam SUM11_mapped.sam
#+   mv HTAdapter88_mapped.sam SUM12_mapped.sam
#+   mv HTAdapter89_mapped.sam SUM13_mapped.sam
#+   mv HTAdapter90_mapped.sam SUM15_mapped.sam
#+   mv HTAdapter91_mapped.sam SUM16_mapped.sam
#+   mv HTAdapter92_mapped.sam SUM17_mapped.sam
#+   mv HTAdapter93_mapped.sam SUM18_mapped.sam
#+   mv HTAdapter94_mapped.sam SUM19_mapped.sam
#+   mv HTAdapter95_mapped.sam SUM20_mapped.sam
#+   mv HTAdapter96_mapped.sam SUM21_mapped.sam
#+   
#+   ## -- ADD READGROUPS -- ##
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_mapped.sam" O="ELF01_IDed.sam" RGPU="ELF" RGSM="ELF_54-38" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_mapped.sam" O="ELF02_IDed.sam" RGPU="ELF" RGSM="ELF_RA607" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_mapped.sam" O="ELF03_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_mapped.sam" O="ELF04_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_mapped.sam" O="ELF05_IDed.sam" RGPU="ELF" RGSM="ELF_558-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_mapped.sam" O="ELF06_IDed.sam" RGPU="ELF" RGSM="ELF_565-13" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_mapped.sam" O="ELF07_IDed.sam" RGPU="ELF" RGSM="ELF_5651-15" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_mapped.sam" O="ELF08_IDed.sam" RGPU="ELF" RGSM="ELF_568-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_mapped.sam" O="ELF09_IDed.sam" RGPU="ELF" RGSM="ELF_546-03" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_mapped.sam" O="ELF10_IDed.sam" RGPU="ELF" RGSM="ELF_547-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_mapped.sam" O="ELF11_IDed.sam" RGPU="ELF" RGSM="ELF_6050-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_mapped.sam" O="ELF12_IDed.sam" RGPU="ELF" RGSM="ELF_622-5" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_mapped.sam" O="ELF13_IDed.sam" RGPU="ELF" RGSM="ELF_636-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_mapped.sam" O="ELF14_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_mapped.sam" O="ELF15_IDed.sam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_mapped.sam" O="ELF16_IDed.sam" RGPU="ELF" RGSM="ELF_RP54" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_mapped.sam" O="MAH01_IDed.sam" RGPU="MAH" RGSM="MAH_35-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_mapped.sam" O="MAH02_IDed.sam" RGPU="MAH" RGSM="MAH_520-02" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_mapped.sam" O="MAH03_IDed.sam" RGPU="MAH" RGSM="MAH_539-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_mapped.sam" O="MAH04_IDed.sam" RGPU="MAH" RGSM="MAH_540-07" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_mapped.sam" O="MAH05_IDed.sam" RGPU="MAH" RGSM="MAH_541-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_mapped.sam" O="MAH06_IDed.sam" RGPU="MAH" RGSM="MAH_542-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_mapped.sam" O="MAH07_IDed.sam" RGPU="MAH" RGSM="MAH_543-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_mapped.sam" O="MAH08_IDed.sam" RGPU="MAH" RGSM="MAH_RA620" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_mapped.sam" O="MAH09_IDed.sam" RGPU="MAH" RGSM="MAH_RA641" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_mapped.sam" O="MAH10_IDed.sam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_mapped.sam" O="MAH11_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_mapped.sam" O="MAH12_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_mapped.sam" O="MAH13_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_mapped.sam" O="MAH14_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_mapped.sam" O="MAH15_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_mapped.sam" O="MAR01_IDed.sam" RGPU="MAR" RGSM="MAR_RP6771" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_mapped.sam" O="MAR02_IDed.sam" RGPU="MAR" RGSM="MAR_RP686" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_mapped.sam" O="MAR03_IDed.sam" RGPU="MAR" RGSM="MAR_368- 2" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_mapped.sam" O="MAR04_IDed.sam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_mapped.sam" O="MAR05_IDed.sam" RGPU="MAR" RGSM="MAR_RA379" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_mapped.sam" O="MAR06_IDed.sam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_mapped.sam" O="MAR07_IDed.sam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_mapped.sam" O="MAR08_IDed.sam" RGPU="MAR" RGSM="MAR_RA605" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_mapped.sam" O="MAR09_IDed.sam" RGPU="MAR" RGSM="MAR_RA630" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_mapped.sam" O="MAR10_IDed.sam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_mapped.sam" O="PAP01_IDed.sam" RGPU="PAP" RGSM="PAP_561-4" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_mapped.sam" O="PAP02_IDed.sam" RGPU="PAP" RGSM="PAP_562-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_mapped.sam" O="PAP03_IDed.sam" RGPU="PAP" RGSM="PAP_564-2" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_mapped.sam" O="PAP04_IDed.sam" RGPU="PAP" RGSM="PAP_569-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_mapped.sam" O="PAP05_IDed.sam" RGPU="PAP" RGSM="PAP_5010-05" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_mapped.sam" O="PAP06_IDed.sam" RGPU="PAP" RGSM="PAP_5021-02" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_mapped.sam" O="PAP07_IDed.sam" RGPU="PAP" RGSM="PAP_505-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_mapped.sam" O="PAP08_IDed.sam" RGPU="PAP" RGSM="PAP_516-02" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_mapped.sam" O="PAP09_IDed.sam" RGPU="PAP" RGSM="PAP_5241-01" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_mapped.sam" O="PAP10_IDed.sam" RGPU="PAP" RGSM="PAP_28-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_mapped.sam" O="PAP11_IDed.sam" RGPU="PAP" RGSM="PAP_RA434" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_mapped.sam" O="PAP12_IDed.sam" RGPU="PAP" RGSM="PAP_RA647" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_mapped.sam" O="PAP13_IDed.sam" RGPU="PAP" RGSM="PAP_RA648" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_mapped.sam" O="PAP14_IDed.sam" RGPU="PAP" RGSM="PAP_RA649" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_mapped.sam" O="PAP15_IDed.sam" RGPU="PAP" RGSM="PAP_RP16" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_mapped.sam" O="STO01_IDed.sam" RGPU="STO" RGSM="STO_RP697" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_mapped.sam" O="STO02_IDed.sam" RGPU="STO" RGSM="STO_RP698" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_mapped.sam" O="STO03_IDed.sam" RGPU="STO" RGSM="STO_RA530" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_mapped.sam" O="STO04_IDed.sam" RGPU="STO" RGSM="STO_RA549" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_mapped.sam" O="STO05_IDed.sam" RGPU="STO" RGSM="STO_RP22" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_mapped.sam" O="STO06_IDed.sam" RGPU="STO" RGSM="STO_271" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_mapped.sam" O="STO07_IDed.sam" RGPU="STO" RGSM="STO_275" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_mapped.sam" O="STO08_IDed.sam" RGPU="STO" RGSM="STO_43-113" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_mapped.sam" O="STO09_IDed.sam" RGPU="STO" RGSM="STO_43-115" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_mapped.sam" O="STO10_IDed.sam" RGPU="STO" RGSM="STO_43-125" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_mapped.sam" O="STO11_IDed.sam" RGPU="STO" RGSM="STO_43-76" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_mapped.sam" O="STO12_IDed.sam" RGPU="STO" RGSM="STO_43-77" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_mapped.sam" O="STO13_IDed.sam" RGPU="STO" RGSM="STO_43-78" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_mapped.sam" O="STO14_IDed.sam" RGPU="STO" RGSM="STO_43-74" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_mapped.sam" O="STO15_IDed.sam" RGPU="STO" RGSM="STO_43-75" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_mapped.sam" O="STO17_IDed.sam" RGPU="STO" RGSM="STO_33-10" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_mapped.sam" O="STO18_IDed.sam" RGPU="STO" RGSM="STO_560-2" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_mapped.sam" O="STO19_IDed.sam" RGPU="STO" RGSM="STO_528-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_mapped.sam" O="STO20_IDed.sam" RGPU="STO" RGSM="STO_559-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_mapped.sam" O="STO21_IDed.sam" RGPU="STO" RGSM="STO_33-272" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_mapped.sam" O="SUM01_IDed.sam" RGPU="SUM" RGSM="SUM_41-103" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_mapped.sam" O="SUM02_IDed.sam" RGPU="SUM" RGSM="SUM_41-104" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_mapped.sam" O="SUM03_IDed.sam" RGPU="SUM" RGSM="SUM_RA645" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_mapped.sam" O="SUM04_IDed.sam" RGPU="SUM" RGSM="SUM_RP31" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_mapped.sam" O="SUM05_IDed.sam" RGPU="SUM" RGSM="SUM_53-255" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_mapped.sam" O="SUM06_IDed.sam" RGPU="SUM" RGSM="SUM_53-252" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_mapped.sam" O="SUM07_IDed.sam" RGPU="SUM" RGSM="SUM_53-253" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_mapped.sam" O="SUM08_IDed.sam" RGPU="SUM" RGSM="SUM_53-256" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_mapped.sam" O="SUM09_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_mapped.sam" O="SUM10_IDed.sam" RGPU="SUM" RGSM="SUM_644-2" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_mapped.sam" O="SUM11_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_mapped.sam" O="SUM12_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_mapped.sam" O="SUM13_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_mapped.sam" O="SUM15_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_mapped.sam" O="SUM16_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_mapped.sam" O="SUM17_IDed.sam" RGPU="SUM" RGSM="SUM_624-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_mapped.sam" O="SUM18_IDed.sam" RGPU="SUM" RGSM="SUM_625-1" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_mapped.sam" O="SUM19_IDed.sam" RGPU="SUM" RGSM="SUM_626-3" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_mapped.sam" O="SUM20_IDed.sam" RGPU="SUM" RGSM="SUM_643-2" RGPL="illumina" RGLB="SeqCap2012"
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_mapped.sam" O="SUM21_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL="illumina" RGLB="SeqCap2012"
#+   
#+   # RNA Seq
#+   cd $WorkingDirectory/mappedReadsRNA
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629651_mapped.sam" O="SRR629651_IDed.sam" RGPU="SUM" RGSM="SAMN01823448" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629599_mapped.sam" O="SRR629599_IDed.sam" RGPU="SUM" RGSM="SAMN01823447" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629652_mapped.sam" O="SRR629652_IDed.sam" RGPU="SUM" RGSM="SAMN01823449" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629653_mapped.sam" O="SRR629653_IDed.sam" RGPU="SUM" RGSM="SAMN01823450" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629654_mapped.sam" O="SRR629654_IDed.sam" RGPU="SUM" RGSM="SAMN01823451" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629655_mapped.sam" O="SRR629655_IDed.sam" RGPU="SUM" RGSM="SAMN01823452" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629656_mapped.sam" O="SRR629656_IDed.sam" RGPU="SUM" RGSM="SAMN01823453" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629657_mapped.sam" O="SRR629657_IDed.sam" RGPU="SUM" RGSM="SAMN01823454" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629658_mapped.sam" O="SRR629658_IDed.sam" RGPU="SUM" RGSM="SAMN01823455" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629659_mapped.sam" O="SRR629659_IDed.sam" RGPU="SUM" RGSM="SAMN01823456" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629661_mapped.sam" O="SRR629661_IDed.sam" RGPU="SUM" RGSM="SAMN01823457" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629660_mapped.sam" O="SRR629660_IDed.sam" RGPU="SUM" RGSM="SAMN01823458" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629662_mapped.sam" O="SRR629662_IDed.sam" RGPU="SUM" RGSM="SAMN01823459" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629664_mapped.sam" O="SRR629664_IDed.sam" RGPU="SUM" RGSM="SAMN01823460" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629663_mapped.sam" O="SRR629663_IDed.sam" RGPU="SUM" RGSM="SAMN01823461" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629665_mapped.sam" O="SRR629665_IDed.sam" RGPU="SUM" RGSM="SAMN01823462" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629666_mapped.sam" O="SRR629666_IDed.sam" RGPU="SUM" RGSM="SAMN01823463" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629667_mapped.sam" O="SRR629667_IDed.sam" RGPU="SUM" RGSM="SAMN01823464" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629668_mapped.sam" O="SRR629668_IDed.sam" RGPU="SUM" RGSM="SAMN01823465" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629669_mapped.sam" O="SRR629669_IDed.sam" RGPU="SUM" RGSM="SAMN01823466" RGPL="illumina" RGLB="RNASeq2012" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497737_mapped.sam" O="SRR497737_IDed.sam" RGPU="SUM" RGSM="SAMN00996377" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497739_mapped.sam" O="SRR497739_IDed.sam" RGPU="SUM" RGSM="SAMN00996379" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497740_mapped.sam" O="SRR497740_IDed.sam" RGPU="SUM" RGSM="SAMN00996380" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497741_mapped.sam" O="SRR497741_IDed.sam" RGPU="SUM" RGSM="SAMN00996381" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497742_mapped.sam" O="SRR497742_IDed.sam" RGPU="SUM" RGSM="SAMN00996382" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497743_mapped.sam" O="SRR497743_IDed.sam" RGPU="SUM" RGSM="SAMN00996383" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497744_mapped.sam" O="SRR497744_IDed.sam" RGPU="SUM" RGSM="SAMN00996384" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497745_mapped.sam" O="SRR497745_IDed.sam" RGPU="SUM" RGSM="SAMN00996385" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497746_mapped.sam" O="SRR497746_IDed.sam" RGPU="SUM" RGSM="SAMN00996386" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497747_mapped.sam" O="SRR497747_IDed.sam" RGPU="SUM" RGSM="SAMN00996387" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497748_mapped.sam" O="SRR497748_IDed.sam" RGPU="SUM" RGSM="SAMN00996388" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497749_mapped.sam" O="SRR497749_IDed.sam" RGPU="SUM" RGSM="SAMN00996389" RGPL="illumina" RGLB="RNASeq2008" 
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497738_mapped.sam" O="SRR497738_IDed.sam" RGPU="SUM" RGSM="SAMN00996378" RGPL="illumina" RGLB="RNASeq2008" 
#+   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   +#

# -- PREPARE FOR SNP CALLING -- ##
# •••••••••••• Functions for obtaining SNP convergence in bqsr ••••••••••••• #
# •••••••••• The input parameter is the current replicate number ••••••••••• #
# ••••••• This approach is called "bootstrapping" by GATK developers ••••••• #
# •••••••••• Even though there is no sampling with replacement... •••••••••• #
function use-HaplotypeCaller {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    -O "$i"_"$1".g.vcf \
    -ERC GVCF
  echo "$i"$'\t'"$i""_""$1"".g.vcf" >> cohort.sample_map_"$1"
done<$WorkingDirectory/mappedReads"$2"/samList
}

function get-just-SNPs {
    # Combine GVCFs
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" GenomicsDBImport \
    	--sample-name-map cohort.sample_map_"$1" \
    	--genomicsdb-workspace-path SNP_database_"$1" \
    	--reader-threads 5 \
    	--intervals genome.intervals
    # Joint genotyping
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" GenotypeGVCFs \
    	-R $WorkingDirectory/References/TelagGenome.fasta \
    	-V gendb://SNP_database_"$1" \
    	-O genotyped_"$1".vcf 
# Get SNPs
/tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" SelectVariants \
	-R $WorkingDirectory/References/TelagGenome.fasta \
	-V genotyped_"$1".vcf \
	--select-type-to-include SNP \
	-O JustSNPs_"$1".vcf
}

function use-BaseRecalibrator {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" BaseRecalibrator \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    --known-sites filtered_0.vcf \
    -O "$i"_recalibration_"$1".table
done<$WorkingDirectory/mappedReads"$2"/samList
}

function use-AnalyzeCovariates {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" AnalyzeCovariates \
    -before "$i"_recalibration_"$1".table \
    -after "$i"_recalibration_"$2".table \
    -plots "$i"_recalComparison_"$1".pdf \
    -csv "$i"_recalComparison_"$1".csv
done<$WorkingDirectory/mappedReads"$3"/samList
}

function use-BQSR {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" ApplyBQSR \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    --bqsr-recal-file "$i"_recalibration_"$1".table \
    -O "$i"_"$2".bam
done<$WorkingDirectory/mappedReads"$3"/samList
}

function combine-VCF {
    #+ cp $WorkingDirectory/GATKDNA/JustSNPs_2.vcf $WorkingDirectory/variantFiltration/JustSNPs_DNA.vcf
    #+ cp $WorkingDirectory/GATKRNA/JustSNPs_2.vcf $WorkingDirectory/variantFiltration/JustSNPs_RNA.vcf
    #+ cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/variantFiltration
    #+ # Change names:
    #+ # DNA
    #+ sed -i.bak "s/SRR497737/ELF52309/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497738/ELF52319/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497739/ELF517112/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497740/PAP50905/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497741/PAP51303/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497742/ELF52503/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497743/ELF54403/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497744/PAP53202/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497745/PAP53307/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497746/PAP53205/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497747/ELF517102/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR497749/PAP50403/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629599/MAH6372/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629651/MAR6271/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629652/MAR6287/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629653/MAR276/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629654/MAR6299/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629655/NAM2064/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629656/MAR6326/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629658/NAM60603/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629659/MAR6463/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629660/NAM6193/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629661/NAM6153/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629662/MAR6099/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629663/NAM6161/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629664/MAR6311/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629665/MAR6341/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629666/MAH252/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629667/MAR6503/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629668/MAH2811/" JustSNPs_RNA.vcf
    #+ sed -i.bak "s/SRR629669/MAH6084/" JustSNPs_RNA.vcf

    #+ # RNA
    #+ sed -i.bak "s/ELF01/ELFRA567/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF02/ELFRA607/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF03/ELFRP540/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF04/ELFRP541/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF05/ELF5581/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF06/ELF56513/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF07/ELF565115/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF08/ELF5683/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF09/ELF54603/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF10/ELF54701/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF11/ELF60503/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF12/ELF6225/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF13/ELF6363/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF14/ELFRA567dup/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF15/ELFRP54.1/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/ELF16/ELFRP54/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH01/MAH351/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH02/MAH52002/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH03/MAH53901/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH04/MAH54007/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH05/MAH54101/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH06/MAH54201/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH07/MAH54301/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH08/MAHRA620/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH09/MAHRA641/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH10/MAHRP47.7/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH11/MAHRP4710/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH12/MAHRP475/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH13/MAHRP478/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH14/MAHRP479/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAH15/MAHRP476/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR01/MARRP6771/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR02/MARRP686/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR03/MAR3682/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR04/MAR3664/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR05/MARRA379/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR06/MAR1794/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR07/MAR1572T/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR08/MARRA605/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR09/MARRA630/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/MAR10/MAR2633/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP01/PAP5614/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP02/PAP5621/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP03/PAP5642/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP04/PAP5693/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP05/PAP501005/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP06/PAP502102/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP07/PAP50501/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP08/PAP51602/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP09/PAP524101/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP10/PAP281/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP11/PAPRA434/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP12/PAPRA647/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP13/PAPRA648/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP14/PAPRA649/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/PAP15/PAPRP16/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO01/STORP697/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO02/STORP698/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO03/STORA530/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO04/STORA549/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO05/STORP22/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO06/STO271/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO07/STO275/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO08/STO43113/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO09/STO43115/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO10/STO43125/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO11/STO4376/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO12/STO4377/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO13/STO4378/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO14/STO4374/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO15/STO4375/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO17/STO3310/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO18/STO5602/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO19/STO5283/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO20/STO5591/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/STO21/STO33272/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM01/SUM41103/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM02/SUM41104/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM03/SUMRA645/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM04/SUMRP31/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM05/SUM53255/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM06/SUM53252/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM07/SUM53253/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM08/SUM53256/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM09/SUMRP531/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM10/SUM6442/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM11/SUMRP5313/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM12/SUMRP5314/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM13/SUMRP535/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM15/SUMRP538/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM16/SUMRP539/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM17/SUM6241/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM18/SUM6251/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM19/SUM6263/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM20/SUM6432/" JustSNPs_DNA.vcf
    #+ sed -i.bak "s/SUM21/SUMRP534/" JustSNPs_DNA.vcf
    #+ cp $WorkingDirectory/variantFiltration/JustSNPs_DNA.vcf $WorkingDirectory/variantFiltration/JustSNPs_DNA_copy.vcf
    #+ cp $WorkingDirectory/variantFiltration/JustSNPs_RNA.vcf $WorkingDirectory/variantFiltration/JustSNPs_RNA_copy.vcf
    cd $WorkingDirectory/variantFiltration/
    #+ bgzip JustSNPs_DNA.vcf
    #+ bgzip JustSNPs_RNA.vcf
    #+ bcftools index JustSNPs_DNA.vcf.gz
    #+ bcftools index JustSNPs_RNA.vcf.gz
    #+ # Find variant intersections between datasets
    #+ bcftools isec JustSNPs_DNA.vcf.gz JustSNPs_RNA.vcf.gz -p isec
    #+ # Move into isec directory
    #+ ## -- This directory contains the output from the isec command
    #+ ## -- Look at the README to understand the files output from isec
    #+ ## -- 0002 is the records from just DNA shared by both DNA and RNA
    #+ ## -- 0003 is the records from just RNA shared by both DNA and RNA
    cd isec
    #+ bgzip 0002.vcf
    #+ bgzip 0003.vcf
    #+ bcftools index 0002.vcf.gz
    #+ bcftools index 0003.vcf.gz
    #+ # Create merged VCF file with the intersected sites
    #+ bcftools merge 0002.vcf.gz 0003.vcf.gz -O v -o Merged.vcf
    #+ vcf2bed < Merged.vcf > isec.bed
    vcftools --gzvcf 0003.vcf.gz --bed isec.bed --out IsecedRNA.vcf --recode --keep-INFO-all
}

function annotateVariants {
  ## ANNOTATE FILTERED VARIANT FILES FOR SEQCAP DATA
  ## ANNOTATE FILTERED VARIANT FILES FOR SEQCAP DATA
  cd $WorkingDirectory/References
  # Make blast database from T. elegans genome
  makeblastdb -in TelagGenome.fasta -parse_seqids -dbtype nucl -out Genome.db
  # Extract IILS genes from exons used for probe design
  python ~/SeqCap/pythonScripts/filterExons.py Exons.fa IILS.txt "$3"TargetGenes.fa log.txt
  # Use Blast with Exons.fa (the exons used for probe design) to filter the genome
  blastn -db Genome.db -query "$3"TargetGenes.fa -outfmt "7 qseqid sseqid evalue qstart qend sstart send" -out BlastResultsIILS.txt
  # Delete "^#" lines from blast output
  cp BlastResultsIILS.txt BlastResultsIILS_original.txt
  sed -i.bak '/^#/d' BlastResultsIILS.txt
  sed -i.bak "s/ref|//" BlastResultsIILS.txt
  sed -i.bak "s/|//" BlastResultsIILS.txt
  # Use filtered genome results (blast output) to pull out targeted genes and create filtered gff
  python ~/SeqCap/pythonScripts/shrinkGFF_v2.py BlastResultsIILS.txt TelagGenome.gff "$3"TargetGenes.gff log.txt
  # Use bedops to convert gff to bed
  gff2bed < "$3"TargetGenes.gff > "$3"TargetGenes.bed
  # bgzip bed file
  bgzip "$3"TargetGenes.bed "$3"TargetGenes.bed.gz
  # tabix index .bed.gz file
  tabix -p bed "$3"TargetGenes.bed.gz
  # Annotate SNP file
  cd $WorkingDirectory/variantFiltration
  bcftools annotate \
  	-a $WorkingDirectory/References/"$3"TargetGenes.bed.gz \
  	-c CHROM,FROM,TO,GENE \
        -o "$2"Init.vcf \
  	-O v \
  	-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  	"$1".vcf
  awk '/^#|GENE=/' "$2"Init.vcf > "$2".vcf
}  

function plotVariants {
source /home/rlk0015/miniconda3/etc/profile.d/conda.sh
conda activate vcfEnv
python ~/SeqCap/pythonScripts/vcf2table.py "$1" "$1"_table
~/SeqCap/RScripts/plotvcftable.R -I "$1"_table -O "$1"_plot.pdf
}

function initial-VariantFiltration {
# Filter variants
	# Provide expression for filtration:
	  # FS = FisherStrand *NOT COMMON IN GATK4.1.7.0
	    # Phred-scaled probability strand bias exists at site 
	  # QD = quality of depth *NOT COMMON IN GATK4.1.7.0
	    # Variant confidence divided by raw depth 
	  # DP = depth of coverage *NOT COMMON IN GATK4.1.7.0
	    # Raw depth. I don't use here, because QD accounts for this 
	  # MQ = RMS mapping quality *NOT COMMON IN GATK4.1.7.0
	    # Root mean square mq over all the reads at the site 
	  # MQRankSum = Mapping quality rank sum test *NOT COMMON IN GATK4.1.7.0
	    # Compares mapping qualities of reads supporting the ref allele to those supporting the alt allele. 
	  # ReadPosRankSum = Rank sum test used to assess site position *NOT COMMON IN GATK4.1.7.0
	    # It compares whether the positions of the reference and alternate alleles are different within the reads 
	  # AB = Allele Balance
	    # Depth of covereage for each allele per sample generalized over all samples
	  # MQ0 = Map Quality 0
	    # Number of reads with map quality 0
	  # SOR = StrandOddsRatio
	    # Estimate strand bias without penalizing reads that occur at the end of exons *FS is biased to penalize
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" VariantFiltration \
	-R $WorkingDirectory/References/TelagGenome.fasta \
	-V "$1".vcf \
        -O "$2"Init.vcf \
	--filter-name "SOR" \
	--filter-expression "SOR > 3.0" \
	--filter-name "QD" \
	--filter-expression "QD < 2.0" \
	--filter-name "MQ" \
	--filter-expression "MQ < 40.0" \
	--filter-name "MQRankSum" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "FS" \
	--filter-expression "FS > 60.0" \
	--filter-name "ReadPosRankSum" \
	--filter-expression "ReadPosRankSum < -5.0"
  awk '/^#/||$7=="PASS"' "$2"Init.vcf > "$2".vcf
} 

function hard-VariantFiltration {
  # Step 1: I changed "VariantFiltration" to filter out SNPs with DP < 20:
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" VariantFiltration -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fasta -V "$1".vcf -O "$2"_HardFilterStep1Init.vcf --filter-name "DP" --filter-expression "DP < 20"
    awk '/^#/||$7=="PASS"' "$2"_HardFilterStep1Init.vcf > "$2"_HardFilterStep1.vcf
  # Step 2: Get rid of multiallelic SNPs (more than 2 alleles):
    bcftools view -m2 -M2 -v snps "$2"_HardFilterStep1.vcf > "$2"_HardFilterStep2.vcf
  # Step 3: Get rid of low-frequency alleles:
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" VariantFiltration -R /scratch/rlk0015/Telag/May2020/WorkingDirectory/References/TelagGenome.fasta -V "$2"_HardFilterStep2.vcf -O "$2"_HardFilterStep3Init.vcf --filter-name "AF" --filter-expression "AF < 0.05"
    awk '/^#/||$7=="PASS"' "$2"_HardFilterStep3Init.vcf > "$2"_HardFilterStep3.vcf
  # Step 4: Get rid of low-quality (mean) genotyping:
    bcftools view  -i  'MIN(FMT/GQ>20)'   "$2"_HardFilterStep3.vcf > "$2"_HardFilterStep4.vcf
  # Step 5: Remove low-depth genotyping:
    vcftools --minDP 20 --vcf "$2"_HardFilterStep4.vcf --recode --recode-INFO-all --out "$2"_HardFilterStep5.vcf
    mv "$2"_HardFilterStep5.vcf.recode.vcf "$2"_HardFilterStep5.vcf
  # Step 6: Remove individual sites with low genotyping
    vcftools --max-missing 0.7 --vcf "$2"_HardFilterStep5.vcf  --recode --recode-INFO-all --out "$2"_HardFilterStep6.vcf
    mv "$2"_HardFilterStep6.vcf.recode.vcf "$2"_HardFiltered.vcf
}
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• #
#+ +++++++++++++++++++++++++++++++ DONE ++++++++++++++++++++++++++++++++ #
#+ # INDEX REF TO USE FOR SNP CALLING
#+ samtools faidx $WorkingDirectory/References/TelagGenome.fasta
#+ java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
#+       R= $WorkingDirectory/References/TelagGenome.fasta \
#+       O= $WorkingDirectory/References/TelagGenome.dict 
#+ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# ============================================================================================ #
# ============================================================================================ #
# =====================    1st Round of SNP calling for SeqCap data    ======================= #
# ============================================================================================ #
# ============================================================================================ #
# == In this section mapped data will be prepped and analyzed, and all the required files   == #
# == for GATK functions will be created.                                                    == #
# ============================================================================================ #

#+   +++++++++++++++++++++++++++++ Done July 4: 1994633 +++++++++++++++++++++++++++++++   +#
#+   cd $WorkingDirectory/mappedReadsDNA
#+   # make .sam files list for Samtools processing
#+   ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
#+   
#+   while read i;
#+   do
#+   cd $WorkingDirectory/mappedReadsDNA
#+   ###  convert .sam to .bam & sort  ###
#+   samtools view -@ 2 -bS $WorkingDirectory/mappedReadsDNA/"$i"*_IDed.sam | samtools sort -@ 2 -o $WorkingDirectory/mappedReadsDNA/"$i"_sorted.bam
#+   ###  remove duplicates  ###
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
#+       INPUT=$WorkingDirectory/mappedReadsDNA/"$i"_sorted.bam \
#+       OUTPUT=$WorkingDirectory/GATKDNA/"$i"_0.bam \
#+       METRICS_FILE=DuplicationMetrics \
#+       CREATE_INDEX=true \
#+       VALIDATION_STRINGENCY=SILENT \
#+       MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#+       ASSUME_SORTED=TRUE \
#+       REMOVE_DUPLICATES=TRUE
#+   #index sorted bam
#+   samtools index $WorkingDirectory/GATKDNA/"$i"_0.bam
#+   	
#+   ###  Calculate Mapping Stats  ###
#+   # tally mapped reads & calcuate the stats
#+   samtools idxstats $WorkingDirectory/GATKDNA/"$i"_0.bam > $WorkingDirectory/StatsDNA/"$i"_counts.txt
#+   samtools flagstat $WorkingDirectory/GATKDNA/"$i"_0.bam > $WorkingDirectory/StatsDNA/"$i"_stats.txt
#+   samtools depth $WorkingDirectory/GATKDNA/"$i"_0.bam > $WorkingDirectory/StatsDNA/"$i"_depth.txt
#+   done<$WorkingDirectory/mappedReadsDNA/samList
#+   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   +#

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

#+   ++++++++++++++++++++++++ Done July 4: 1994303 & 1994412 ++++++++++++++++++++++++   +#
#+   ### --------------------- Continue Variant Calling ----------------------- ###
#+   # ============================================================================================ #
#+   # ============================================================================================ #
#+   # =====================    1st Round of SNP calling for RNAseq data    ======================= #
#+   # ============================================================================================ #
#+   # ============================================================================================ #
#+   # == In this section mapped data will be prepped and analyzed, and all the required files   == #
#+   # == for GATK functions will be created.                                                    == #
#+   # ============================================================================================ #
#+   
#+   cd $WorkingDirectory/mappedReadsRNA
#+   # make .sam files list for Samtools processing
#+   ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
#+   while read i;
#+   do
#+   ###  convert .sam to .bam & sort  ###
#+   samtools view -@ 2 -bS $WorkingDirectory/mappedReadsRNA/"$i"*_IDed.sam | samtools sort -@ 2 -o $WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam
#+   ###  remove duplicates  ###
#+   java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
#+       INPUT=$WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam \
#+       OUTPUT=$WorkingDirectory/GATKRNA/"$i"_0.bam \
#+       METRICS_FILE=DuplicationMetrics \
#+       CREATE_INDEX=true \
#+       VALIDATION_STRINGENCY=SILENT \
#+       MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#+       ASSUME_SORTED=TRUE \
#+       REMOVE_DUPLICATES=TRUE
#+   #index sorted bam
#+   samtools index $WorkingDirectory/GATKRNA/"$i"_0.bam
#+   	
#+   ###  Calculate Mapping Stats  ###
#+   # tally mapped reads & calcuate the stats
#+   samtools idxstats $WorkingDirectory/GATKRNA/"$i"_0.bam > $WorkingDirectory/StatsRNA/"$i"_counts.txt
#+   samtools flagstat $WorkingDirectory/GATKRNA/"$i"_0.bam > $WorkingDirectory/StatsRNA/"$i"_stats.txt
#+   samtools depth $WorkingDirectory/GATKRNA/"$i"_0.bam > $WorkingDirectory/StatsRNA/"$i"_depth.txt
#+   done<$WorkingDirectory/mappedReadsRNA/samList
#+   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   +#

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

#+ +++++++++++++++++++++++++++ Completed July 2020 +++++++++++++++++++++++++++ +#
#+ ### --------------------- Continue Variant Calling ----------------------- ###
#+ # Perform bqsr-'bootstraping', SNP calling, and hard filtration on Seq Cap data
#+ cd $WorkingDirectory/GATKDNA
#+ ## Replicate 1
#+ use-HaplotypeCaller 0 DNA
#+ get-just-SNPs 0 
#+ use-VariantFiltration 0
#+ use-BaseRecalibrator 0 DNA
#+ use-BQSR 0 1 DNA
#+ use-HaplotypeCaller 1 DNA
#+ get-just-SNPs 1
#+ use-BaseRecalibrator 1 DNA
#+ use-AnalyzeCovariates 0 1 DNA
## -- Replicate 2
#+ use-BQSR 1 2 DNA
#+ use-HaplotypeCaller 2 DNA
#+ get-just-SNPs 2
#+ use-BaseRecalibrator 2 DNA
#+ use-AnalyzeCovariates 1 2 DNA
#+ ## -- Merge RNA and DNA data
#+ combine-VCF
#+ ## -- Annotate variants
#+ annotateVariants Merged IILS_Annotated IILS
#+ annotateVariants Merged All_Annotated All
#+ ## -- Examine Unfiltered Variants
cd $WorkingDirectory/variantFiltration
#+ plotVariants IILS_Annotated.vcf
#+ plotVariants SeqCap_Annotated.vcf
#+ plotVariants All_Annotated.vcf
#+ ## -- Initial Filter Variants
#+ initial-VariantFiltration IILS_Annotated IILS_InitialFiltered
#+ initial-VariantFiltration SeqCap_Annotated SeqCap_InitialFiltered
#+ initial-VariantFiltration All_Annotated All_InitialFiltered
#+ ## -- Examine Initial Filtered Variants
#+ plotVariants IILS_InitialFiltered.vcf
#+ plotVariants SeqCap_InitialFiltered.vcf
#+ plotVariants All_InitialFiltered.vcf
#+ ## -- Perform Hard Filtering
#+ hard-VariantFiltration IILS_InitialFiltered IILS
#+ hard-VariantFiltration SeqCap_InitialFiltered SeqCap
#+ hard-VariantFiltration All_InitialFiltered All
#+ ## -- Examine Hard Filtered Variants
#+ plotVariants IILS_HardFiltered.vcf
#+ plotVariants SeqCap_HardFiltered.vcf
#+ plotVariants All_HardFiltered.vcf
#+ cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/variantFiltration/isec
#+ vcftools --vcf Merged.vcf --bed isec.bed --out Iseced.vcf --recode --keep-INFO-all
#+ 
#+ # Perform bqsr-'bootstraping', SNP calling, and hard filtration on RNA-seq data
#+ cd $WorkingDirectory/GATKRNA
#+ ## Replicate 1
#+ use-HaplotypeCaller 0 RNA
#+ get-just-SNPs 0
#+ use-VariantFiltration 0
#+ use-BaseRecalibrator 0 RNA
#+ use-BQSR 0 1 RNA
#+ use-HaplotypeCaller 1 RNA
#+ get-just-SNPs 1
#+ use-BaseRecalibrator 1 RNA
#+ use-AnalyzeCovariates 0 1 RNA
#+ ## -- Replicate 2
#+ use-BQSR 1 2 RNA
#+ use-HaplotypeCaller 2 RNA
#+ get-just-SNPs 2
#+ use-BaseRecalibrator 2 RNA
#+ use-AnalyzeCovariates 1 2 RNA
#+ ## -- Filter Variants
combine-VCF
#+ use-VariantFiltration 2
#+ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+ +++++++++++++++++++++++++++  Completed  +++++++++++++++++++++++++++ +#
#+ # copy gff annotated genome to references
#+ cp /home/rlk0015/SeqCap/code/References/T_elegans_genome/latest_assembly_versions/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.gff.gz $WorkingDirectory/References/TelagGenome.gff.gz
#+ cd $WorkingDirectory/References
#+ gunzip TelagGenome.gff.gz
#+ echo "RLK_report: REFERENCE GFF COPY COMPLETE"
#+ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

