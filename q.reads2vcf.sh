#!/bin/sh

#Give job a name
#PBS -N ThamnophisSeqCapPop2

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

#create working directory ###
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads
echo "RLK_report: directory created: /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome with rawReads and cleanReads sub directories"

### copy raw data to scratch ###
cd /scratch/rlk0015/Telag/rawReads
ls HTAdapter*.fastq.gz | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads'

### quality check ###
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads
ls *.fastq.gz | time parallel -j+0 --eta 'fastqc {}'

### copy over adapter file ###
cp /home/rlk0015/SeqCap/code/References/adapters.fa /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads

### paired-end trimming ###
ls | grep "fastq.gz" | cut -d "_" -f 1,2 | sort | uniq > PE_TrimmList
while read i
do
java -jar /tools/trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads/"$i"*_R1*.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/rawReads/"$i"*_R2*.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"_R1_paired.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"_R1_unpaired.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"_R2_paired.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"_R2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:30 MINLEN:36
done<PE_TrimmList

### quality check ###
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/
ls *paired.fastq.gz | time parallel -j+0 --eta 'fastqc {}'

## -- MAP TO REFERENCE -- ##
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads
### copy the reference to mapped reads directory ###
cp /home/rlk0015/SeqCap/code/References/ReducedTranscriptome.fa /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads
echo "RLK_report: reference copy complete"
# index reference
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads
bwa index -p transcriptome -a is ReducedTranscriptome.fa
echo "RLK_report: reference index complete"

# create list with each paired individual
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads
ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/pairedMapList
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads
# while loop through the names in pairedMapList
while read i
do
### map to ref transcriptome ###
bwa mem -t 4 -M transcriptome /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"*_R1_paired.fastq.gz /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/cleanReads/"$i"*_R2_paired.fastq.gz > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/"$i"_doubles.sam
done<pairedMapList
echo "RLK_report: map complete"

### change names to match population sampling ###
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads
mv HTAdapter1_doubles.sam ELF01_doubles.sam
mv HTAdapter2_doubles.sam ELF02_doubles.sam
mv HTAdapter3_doubles.sam ELF03_doubles.sam
mv HTAdapter4_doubles.sam ELF04_doubles.sam
mv HTAdapter5_doubles.sam ELF05_doubles.sam
mv HTAdapter6_doubles.sam ELF06_doubles.sam
mv HTAdapter7_doubles.sam ELF07_doubles.sam
mv HTAdapter8_doubles.sam ELF08_doubles.sam
mv HTAdapter9_doubles.sam ELF09_doubles.sam
mv HTAdapter10_doubles.sam ELF10_doubles.sam
mv HTAdapter11_doubles.sam ELF11_doubles.sam
mv HTAdapter12_doubles.sam ELF12_doubles.sam
mv HTAdapter13_doubles.sam ELF13_doubles.sam
mv HTAdapter14_doubles.sam ELF14_doubles.sam
mv HTAdapter15_doubles.sam ELF15_doubles.sam
mv HTAdapter16_doubles.sam ELF16_doubles.sam
mv HTAdapter17_doubles.sam MAH01_doubles.sam
mv HTAdapter18_doubles.sam MAH02_doubles.sam
mv HTAdapter19_doubles.sam MAH03_doubles.sam
mv HTAdapter20_doubles.sam MAH04_doubles.sam
mv HTAdapter21_doubles.sam MAH05_doubles.sam
mv HTAdapter22_doubles.sam MAH06_doubles.sam
mv HTAdapter23_doubles.sam MAH07_doubles.sam
mv HTAdapter24_doubles.sam MAH08_doubles.sam
mv HTAdapter25_doubles.sam MAH09_doubles.sam
mv HTAdapter26_doubles.sam MAH10_doubles.sam
mv HTAdapter27_doubles.sam MAH11_doubles.sam
mv HTAdapter28_doubles.sam MAH12_doubles.sam
mv HTAdapter29_doubles.sam MAH13_doubles.sam
mv HTAdapter30_doubles.sam MAH14_doubles.sam
mv HTAdapter31_doubles.sam MAH15_doubles.sam
mv HTAdapter32_doubles.sam MAR01_doubles.sam
mv HTAdapter33_doubles.sam MAR02_doubles.sam
mv HTAdapter34_doubles.sam MAR03_doubles.sam
mv HTAdapter35_doubles.sam MAR04_doubles.sam
mv HTAdapter36_doubles.sam MAR05_doubles.sam
mv HTAdapter37_doubles.sam MAR06_doubles.sam
mv HTAdapter38_doubles.sam MAR07_doubles.sam
mv HTAdapter39_doubles.sam MAR08_doubles.sam
mv HTAdapter40_doubles.sam MAR09_doubles.sam
mv HTAdapter41_doubles.sam MAR10_doubles.sam
mv HTAdapter42_doubles.sam PAP01_doubles.sam
mv HTAdapter43_doubles.sam PAP02_doubles.sam
mv HTAdapter44_doubles.sam PAP03_doubles.sam
mv HTAdapter45_doubles.sam PAP04_doubles.sam
mv HTAdapter46_doubles.sam PAP05_doubles.sam
mv HTAdapter47_doubles.sam PAP06_doubles.sam
mv HTAdapter48_doubles.sam PAP07_doubles.sam
mv HTAdapter49_doubles.sam PAP08_doubles.sam
mv HTAdapter50_doubles.sam PAP09_doubles.sam
mv HTAdapter51_doubles.sam PAP10_doubles.sam
mv HTAdapter52_doubles.sam PAP11_doubles.sam
mv HTAdapter53_doubles.sam PAP12_doubles.sam
mv HTAdapter54_doubles.sam PAP13_doubles.sam
mv HTAdapter55_doubles.sam PAP14_doubles.sam
mv HTAdapter56_doubles.sam PAP15_doubles.sam
mv HTAdapter57_doubles.sam STO01_doubles.sam
mv HTAdapter58_doubles.sam STO02_doubles.sam
mv HTAdapter59_doubles.sam STO03_doubles.sam
mv HTAdapter60_doubles.sam STO04_doubles.sam
mv HTAdapter61_doubles.sam STO05_doubles.sam
mv HTAdapter62_doubles.sam STO06_doubles.sam
mv HTAdapter63_doubles.sam STO07_doubles.sam
mv HTAdapter64_doubles.sam STO08_doubles.sam
mv HTAdapter65_doubles.sam STO09_doubles.sam
mv HTAdapter66_doubles.sam STO10_doubles.sam
mv HTAdapter67_doubles.sam STO11_doubles.sam
mv HTAdapter68_doubles.sam STO12_doubles.sam
mv HTAdapter69_doubles.sam STO13_doubles.sam
mv HTAdapter70_doubles.sam STO14_doubles.sam
mv HTAdapter71_doubles.sam STO15_doubles.sam
mv HTAdapter72_doubles.sam STO17_doubles.sam
mv HTAdapter73_doubles.sam STO18_doubles.sam
mv HTAdapter74_doubles.sam STO19_doubles.sam
mv HTAdapter75_doubles.sam STO20_doubles.sam
mv HTAdapter76_doubles.sam STO21_doubles.sam
mv HTAdapter77_doubles.sam SUM01_doubles.sam
mv HTAdapter78_doubles.sam SUM02_doubles.sam
mv HTAdapter79_doubles.sam SUM03_doubles.sam
mv HTAdapter80_doubles.sam SUM04_doubles.sam
mv HTAdapter81_doubles.sam SUM05_doubles.sam
mv HTAdapter82_doubles.sam SUM06_doubles.sam
mv HTAdapter83_doubles.sam SUM07_doubles.sam
mv HTAdapter84_doubles.sam SUM08_doubles.sam
mv HTAdapter85_doubles.sam SUM09_doubles.sam
mv HTAdapter86_doubles.sam SUM10_doubles.sam
mv HTAdapter87_doubles.sam SUM11_doubles.sam
mv HTAdapter88_doubles.sam SUM12_doubles.sam
mv HTAdapter89_doubles.sam SUM13_doubles.sam
mv HTAdapter90_doubles.sam SUM15_doubles.sam
mv HTAdapter91_doubles.sam SUM16_doubles.sam
mv HTAdapter92_doubles.sam SUM17_doubles.sam
mv HTAdapter93_doubles.sam SUM18_doubles.sam
mv HTAdapter94_doubles.sam SUM19_doubles.sam
mv HTAdapter95_doubles.sam SUM20_doubles.sam
mv HTAdapter96_doubles.sam SUM21_doubles.sam

## -- ADD READGROUPS -- ##
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_doubles.sam" O="ELF01_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_doubles.sam" O="ELF02_IDed.sam" RGPU="ELF" RGSM="ELF_RA607" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_doubles.sam" O="ELF03_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_doubles.sam" O="ELF04_IDed.sam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_doubles.sam" O="ELF05_IDed.sam" RGPU="ELF" RGSM="ELF_558-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_doubles.sam" O="ELF06_IDed.sam" RGPU="ELF" RGSM="ELF_565-13" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_doubles.sam" O="ELF07_IDed.sam" RGPU="ELF" RGSM="ELF_5651-15" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_doubles.sam" O="ELF08_IDed.sam" RGPU="ELF" RGSM="ELF_568-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_doubles.sam" O="ELF09_IDed.sam" RGPU="ELF" RGSM="ELF_546-03" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_doubles.sam" O="ELF10_IDed.sam" RGPU="ELF" RGSM="ELF_547-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_doubles.sam" O="ELF11_IDed.sam" RGPU="ELF" RGSM="ELF_6050-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_doubles.sam" O="ELF12_IDed.sam" RGPU="ELF" RGSM="ELF_622-5" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_doubles.sam" O="ELF13_IDed.sam" RGPU="ELF" RGSM="ELF_636-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_doubles.sam" O="ELF14_IDed.sam" RGPU="ELF" RGSM="ELF_RA567" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_doubles.sam" O="ELF15_IDed.sam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_doubles.sam" O="ELF16_IDed.sam" RGPU="ELF" RGSM="ELF_RP54" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_doubles.sam" O="MAH01_IDed.sam" RGPU="MAH" RGSM="MAH_35-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_doubles.sam" O="MAH02_IDed.sam" RGPU="MAH" RGSM="MAH_520-02" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_doubles.sam" O="MAH03_IDed.sam" RGPU="MAH" RGSM="MAH_539-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_doubles.sam" O="MAH04_IDed.sam" RGPU="MAH" RGSM="MAH_540-07" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_doubles.sam" O="MAH05_IDed.sam" RGPU="MAH" RGSM="MAH_541-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_doubles.sam" O="MAH06_IDed.sam" RGPU="MAH" RGSM="MAH_542-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_doubles.sam" O="MAH07_IDed.sam" RGPU="MAH" RGSM="MAH_543-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_doubles.sam" O="MAH08_IDed.sam" RGPU="MAH" RGSM="MAH_RA620" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_doubles.sam" O="MAH09_IDed.sam" RGPU="MAH" RGSM="MAH_RA641" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_doubles.sam" O="MAH10_IDed.sam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_doubles.sam" O="MAH11_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_doubles.sam" O="MAH12_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_doubles.sam" O="MAH13_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_doubles.sam" O="MAH14_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_doubles.sam" O="MAH15_IDed.sam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_doubles.sam" O="MAR01_IDed.sam" RGPU="MAR" RGSM="MAR_RP6771" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_doubles.sam" O="MAR02_IDed.sam" RGPU="MAR" RGSM="MAR_RP686" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_doubles.sam" O="MAR03_IDed.sam" RGPU="MAR" RGSM="MAR_368- 2" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_doubles.sam" O="MAR04_IDed.sam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_doubles.sam" O="MAR05_IDed.sam" RGPU="MAR" RGSM="MAR_RA379" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_doubles.sam" O="MAR06_IDed.sam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_doubles.sam" O="MAR07_IDed.sam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_doubles.sam" O="MAR08_IDed.sam" RGPU="MAR" RGSM="MAR_RA605" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_doubles.sam" O="MAR09_IDed.sam" RGPU="MAR" RGSM="MAR_RA630" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_doubles.sam" O="MAR10_IDed.sam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_doubles.sam" O="PAP01_IDed.sam" RGPU="PAP" RGSM="PAP_561-4" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_doubles.sam" O="PAP02_IDed.sam" RGPU="PAP" RGSM="PAP_562-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_doubles.sam" O="PAP03_IDed.sam" RGPU="PAP" RGSM="PAP_564-2" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_doubles.sam" O="PAP04_IDed.sam" RGPU="PAP" RGSM="PAP_569-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_doubles.sam" O="PAP05_IDed.sam" RGPU="PAP" RGSM="PAP_5010-05" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_doubles.sam" O="PAP06_IDed.sam" RGPU="PAP" RGSM="PAP_5021-02" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_doubles.sam" O="PAP07_IDed.sam" RGPU="PAP" RGSM="PAP_505-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_doubles.sam" O="PAP08_IDed.sam" RGPU="PAP" RGSM="PAP_516-02" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_doubles.sam" O="PAP09_IDed.sam" RGPU="PAP" RGSM="PAP_5241-01" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_doubles.sam" O="PAP10_IDed.sam" RGPU="PAP" RGSM="PAP_28-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_doubles.sam" O="PAP11_IDed.sam" RGPU="PAP" RGSM="PAP_RA434" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_doubles.sam" O="PAP12_IDed.sam" RGPU="PAP" RGSM="PAP_RA647" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_doubles.sam" O="PAP13_IDed.sam" RGPU="PAP" RGSM="PAP_RA648" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_doubles.sam" O="PAP14_IDed.sam" RGPU="PAP" RGSM="PAP_RA649" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_doubles.sam" O="PAP15_IDed.sam" RGPU="PAP" RGSM="PAP_RP16" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_doubles.sam" O="STO01_IDed.sam" RGPU="STO" RGSM="STO_RP697" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_doubles.sam" O="STO02_IDed.sam" RGPU="STO" RGSM="STO_RP698" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_doubles.sam" O="STO03_IDed.sam" RGPU="STO" RGSM="STO_RA530" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_doubles.sam" O="STO04_IDed.sam" RGPU="STO" RGSM="STO_RA549" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_doubles.sam" O="STO05_IDed.sam" RGPU="STO" RGSM="STO_RP22" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_doubles.sam" O="STO06_IDed.sam" RGPU="STO" RGSM="STO_RP33" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_doubles.sam" O="STO07_IDed.sam" RGPU="STO" RGSM="STO_RP33" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_doubles.sam" O="STO08_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_doubles.sam" O="STO09_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_doubles.sam" O="STO10_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_doubles.sam" O="STO11_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_doubles.sam" O="STO12_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_doubles.sam" O="STO13_IDed.sam" RGPU="STO" RGSM="STO_NC43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_doubles.sam" O="STO14_IDed.sam" RGPU="STO" RGSM="STO_RP43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_doubles.sam" O="STO15_IDed.sam" RGPU="STO" RGSM="STO_RP43" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_doubles.sam" O="STO17_IDed.sam" RGPU="STO" RGSM="STO_33-10" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_doubles.sam" O="STO18_IDed.sam" RGPU="STO" RGSM="STO_560-2" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_doubles.sam" O="STO19_IDed.sam" RGPU="STO" RGSM="STO_528-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_doubles.sam" O="STO20_IDed.sam" RGPU="STO" RGSM="STO_559-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_doubles.sam" O="STO21_IDed.sam" RGPU="STO" RGSM="STO_RP33" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_doubles.sam" O="SUM01_IDed.sam" RGPU="SUM" RGSM="SUM_RP41" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_doubles.sam" O="SUM02_IDed.sam" RGPU="SUM" RGSM="SUM_RP41" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_doubles.sam" O="SUM03_IDed.sam" RGPU="SUM" RGSM="SUM_RA645" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_doubles.sam" O="SUM04_IDed.sam" RGPU="SUM" RGSM="SUM_RP31" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_doubles.sam" O="SUM05_IDed.sam" RGPU="SUM" RGSM="SUM_RP41" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_doubles.sam" O="SUM06_IDed.sam" RGPU="SUM" RGSM="SUM_RP53" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_doubles.sam" O="SUM07_IDed.sam" RGPU="SUM" RGSM="SUM_RP53" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_doubles.sam" O="SUM08_IDed.sam" RGPU="SUM" RGSM="SUM_RP53" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_doubles.sam" O="SUM09_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_doubles.sam" O="SUM10_IDed.sam" RGPU="SUM" RGSM="SUM_644-2" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_doubles.sam" O="SUM11_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_doubles.sam" O="SUM12_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_doubles.sam" O="SUM13_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_doubles.sam" O="SUM15_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_doubles.sam" O="SUM16_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_doubles.sam" O="SUM17_IDed.sam" RGPU="SUM" RGSM="SUM_624-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_doubles.sam" O="SUM18_IDed.sam" RGPU="SUM" RGSM="SUM_625-1" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_doubles.sam" O="SUM19_IDed.sam" RGPU="SUM" RGSM="SUM_626-3" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_doubles.sam" O="SUM20_IDed.sam" RGPU="SUM" RGSM="SUM_643-2" RGPL:"illumina" RGLB:"SeqCap2012"
java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_doubles.sam" O="SUM21_IDed.sam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL:"illumina" RGLB:"SeqCap2012"
echo "RLK_report: readGrouping complete"

## -- PREPARE FOR SNP CALLING -- ##
# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# make variant directory
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/SNPTables
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/sequenceTables
# copy reference for SNP calling
cp /home/rlk0015/SeqCap/code/References/ReducedTranscriptome.fa /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa
#index ref
samtools faidx /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa
java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
    R= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    O= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/Transcriptome.dict
## ** This didn't work ** ##
# java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
#     R= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
#     O= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/TranscriptomeForGATK.dict
## ********************** ##
# cp /home/rlk0015/SeqCap/code/References/ReducedTranscriptome.fa # /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa
# samtools faidx /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa
# java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
#     R= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
#     O= /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/Transcriptome.dict
# loop through .sam files
while read i;
do
### convert .sam to .bam & sort. ex: HTAdapter84_CACCTTAC_L003_R1.sam –> HTAdapter84_CACCTTAC_L003_R1_sorted.bam ###
samtools view -@ 2 -bS /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/"$i"_IDed.sam | samtools sort -@ 2 -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/"$i"_sorted.bam
### remove duplicates ###
java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
    INPUT=/scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/"$i"_sorted.bam \
    OUTPUT=/scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_dupsRemoved.bam \
    METRICS_FILE=DuplicationMetrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    ASSUME_SORTED=TRUE \
    REMOVE_DUPLICATES=TRUE
#index sorted bam
samtools index "$i"_dupsRemoved.bam
done<samList
### move to GATK directory ###
cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/
### merge .bam files ###
samtools merge /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/dupsRemoved.bam /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/*_dupsRemoved.bam
### index the merged .bam ###
samtools index /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/dupsRemoved.bam
# call indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/dupsRemoved.bam \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsCalled.intervals
# realign indels
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/dupsRemoved.bam \
    -targetIntervals /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsCalled.intervals \
    -LOD 3.0 \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsRealigned.bam
## -- CALL SNPS -- ##
java -Xmx16g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsRealigned.bam \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/calledSNPs.vcf \
    -gt_mode DISCOVERY \
    -ploidy 2 \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# annotate variants
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsRealigned.bam \
    -G StandardAnnotation \
    -V:variant,VCF /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/calledSNPs.vcf \
    -XA SnpEff \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/annotatedSNPs.vcf \
    -rf BadCigar
# annotate indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsRealigned.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/annotatedIndels.vcf \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -rf BadCigar
# mask indels
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -V /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/calledSNPs.vcf \
    --mask /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/annotatedIndels.vcf \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/SNPsMaskedIndels.vcf/ \
    -rf BadCigar
# restrict to high-quality variant calls
cat /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/SNPsMaskedIndels.vcf | grep 'PASS\|^#' > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/qualitySNPs.vcf
# read-backed phasing
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T ReadBackedPhasing \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -I /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/indelsRealigned.bam \
    --variant /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/qualitySNPs.vcf \
    -L /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/qualitySNPs.vcf \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/phasedSNPs.vcf \
    --phaseQualityThresh 20.0 \
    -rf BadCigar
while read i;
do
# VCF for each sample
java -Xmx2g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    --variant /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/phasedSNPs.vcf \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_phasedSNPs.vcf \
    -sn "$i" \
    -rf BadCigar
# make SNPs table
java -Xmx8g -jar /tools/gatk-3.6/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    -V /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_phasedSNPs.vcf \
    -F CHROM -F POS -F QUAL -GF GT -GF DP -GF HP -GF AD \
    -o /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/SNPTables/"$i"_tableSNPs.txt \
    -rf BadCigar
# Add phased SNPs to reference and filter
python /home/rlk0015/SeqCap/seqcap_pop/bin/add_phased_snps_to_seqs_filter.py \
    /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/ReducedTranscriptome.fa \
    /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/SNPTables/"$i"_tableSNPs.txt \
    /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/sequenceTables/"$i"_tableSequences.txt \
    1
# make stats folder
mkdir -p /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/stats
# tally mapped reads & calcuate the stats
samtools idxstats /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/stats/"$i"_counts.txt
samtools flagstat /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/stats/"$i"_stats.txt
samtools depth /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/GATK/"$i"_dupsRemoved.bam > /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/stats/"$i"_depth.txt
done</scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/mappedReads/samList

cd /scratch/rlk0015/Telag/Aug2018/WorkingDirectory_ReducedTranscriptome/stats
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList

# Generate file with all depth information
cat *_depth.txt > Transcriptome_depth.txt
for f in Transcriptome_depth.txt
do
sed -i "s/$/\tTranscriptome/" $f; done

# Create file with avg depth per exon using Randy's python script
#!!!!!!!! CHANGE THIS TO "PER CONTIG" FOR TRANSCRIPTOME (WILL HAVE TO CHANGE CODE)
python /home/rlk0015/SeqCap/pythonScripts/avgDepth.py Transcriptome_depth.txt Transcriptome_avgDepth.txt

# make results directory & move results
mkdir -p /home/rlk0015/SeqCap/Aug2018/stats
mkdir -p /home/rlk0015/SeqCap/Aug2018/counts
mkdir -p /home/rlk0015/SeqCap/Aug2018/avgDepth

# cp results to respective directories
cp *stats.txt /home/rlk0015/SeqCap/Aug2018/stats
cp *counts.txt /home/rlk0015/SeqCap/Aug2018/counts
cp *depth.txt /home/rlk0015/SeqCap/Aug2018/avgDepth

