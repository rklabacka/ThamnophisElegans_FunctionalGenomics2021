#!/bin/sh

#Give job a name
#PBS -N completeExons_mapToExons_t3

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=3:ppn=15,walltime=40:00:00 
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

# #create working directory
# mkdir -p /scratch/rlk0015/Telag/Results_CompleteExonList/MapToExons/WorkingDirectory
# echo "RLK_report: directory created: /scratch/rlk0015/Telag/Results_CompleteExonList/MapToExons/WorkingDirectory"
# 
# # copy over clean SE fastq files
# cd /scratch/rlk0015/Telag/CleanedReads/SingleEnd/
# ls *_trimmed.fastq | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Results_CompleteExonList/MapToExons/WorkingDirectory'
# echo "RLK_report: clean SE fastq files copied to /scratch/rlk0015/Telag/Results_CompleteExonList/MapToExons/WorkingDirectory"
# 
# move to working directory
cd /scratch/rlk0015/Telag/Results_CompleteExonList/MapToExons/WorkingDirectory
# 
# # create list for each trimmed individual
# ls | grep "_trimmed.fastq" |cut -d "_" -f 1,2,3 | sort | uniq > mergeList
# # while loop through the names in mergeList
# while read i
# do
# # merge R1 & R2 files into single file
# cat "$i"_R*_trimmed.fastq > "$i"_merged.fastq
# done<mergeList
# echo "RLK_report: merge complete"

# copy the reference to your current directory
cp /home/rlk0015/SeqCap/code/References/CompleteReferenceList/CompleteExons_DupsRemoved_Randy.fa .
echo "RLK_report: reference copy complete"

# index reference
bwa index -p exons -a is CompleteExons_DupsRemoved_Randy.fa
echo "RLK_report: reference index complete"

# create list with each merged individual
ls | grep "_merged.fastq" |cut -d "_" -f 1 | sort | uniq > mapList
# while loop through the names in mapList
while read i
do
# map to ref
bwa mem -t 4 -M exons "$i"*_merged.fastq > "$i"_singles.sam
done<mapList
echo "RLK_report: map complete"

# map: Tonia's method
mv HTAdapter1_singles.sam ELF01_singles.sam
mv HTAdapter2_singles.sam ELF02_singles.sam
mv HTAdapter3_singles.sam ELF03_singles.sam
mv HTAdapter4_singles.sam ELF04_singles.sam
mv HTAdapter5_singles.sam ELF05_singles.sam
mv HTAdapter6_singles.sam ELF06_singles.sam
mv HTAdapter7_singles.sam ELF07_singles.sam
mv HTAdapter8_singles.sam ELF08_singles.sam
mv HTAdapter9_singles.sam ELF09_singles.sam
mv HTAdapter10_singles.sam ELF10_singles.sam
mv HTAdapter11_singles.sam ELF11_singles.sam
mv HTAdapter12_singles.sam ELF12_singles.sam
mv HTAdapter13_singles.sam ELF13_singles.sam
mv HTAdapter14_singles.sam ELF14_singles.sam
mv HTAdapter15_singles.sam ELF15_singles.sam
mv HTAdapter16_singles.sam ELF16_singles.sam
mv HTAdapter17_singles.sam MAH01_singles.sam
mv HTAdapter18_singles.sam MAH02_singles.sam
mv HTAdapter19_singles.sam MAH03_singles.sam
mv HTAdapter20_singles.sam MAH04_singles.sam
mv HTAdapter21_singles.sam MAH05_singles.sam
mv HTAdapter22_singles.sam MAH06_singles.sam
mv HTAdapter23_singles.sam MAH07_singles.sam
mv HTAdapter24_singles.sam MAH08_singles.sam
mv HTAdapter25_singles.sam MAH09_singles.sam
mv HTAdapter26_singles.sam MAH10_singles.sam
mv HTAdapter27_singles.sam MAH11_singles.sam
mv HTAdapter28_singles.sam MAH12_singles.sam
mv HTAdapter29_singles.sam MAH13_singles.sam
mv HTAdapter30_singles.sam MAH14_singles.sam
mv HTAdapter31_singles.sam MAH15_singles.sam
mv HTAdapter32_singles.sam MAR01_singles.sam
mv HTAdapter33_singles.sam MAR02_singles.sam
mv HTAdapter34_singles.sam MAR03_singles.sam
mv HTAdapter35_singles.sam MAR04_singles.sam
mv HTAdapter36_singles.sam MAR05_singles.sam
mv HTAdapter37_singles.sam MAR06_singles.sam
mv HTAdapter38_singles.sam MAR07_singles.sam
mv HTAdapter39_singles.sam MAR08_singles.sam
mv HTAdapter40_singles.sam MAR09_singles.sam
mv HTAdapter41_singles.sam MAR10_singles.sam
mv HTAdapter42_singles.sam PAP01_singles.sam
mv HTAdapter43_singles.sam PAP02_singles.sam
mv HTAdapter44_singles.sam PAP03_singles.sam
mv HTAdapter45_singles.sam PAP04_singles.sam
mv HTAdapter46_singles.sam PAP05_singles.sam
mv HTAdapter47_singles.sam PAP06_singles.sam
mv HTAdapter48_singles.sam PAP07_singles.sam
mv HTAdapter49_singles.sam PAP08_singles.sam
mv HTAdapter50_singles.sam PAP09_singles.sam
mv HTAdapter51_singles.sam PAP10_singles.sam
mv HTAdapter52_singles.sam PAP11_singles.sam
mv HTAdapter53_singles.sam PAP12_singles.sam
mv HTAdapter54_singles.sam PAP13_singles.sam
mv HTAdapter55_singles.sam PAP14_singles.sam
mv HTAdapter56_singles.sam PAP15_singles.sam
mv HTAdapter57_singles.sam STO01_singles.sam
mv HTAdapter58_singles.sam STO02_singles.sam
mv HTAdapter59_singles.sam STO03_singles.sam
mv HTAdapter60_singles.sam STO04_singles.sam
mv HTAdapter61_singles.sam STO05_singles.sam
mv HTAdapter62_singles.sam STO06_singles.sam
mv HTAdapter63_singles.sam STO07_singles.sam
mv HTAdapter64_singles.sam STO08_singles.sam
mv HTAdapter65_singles.sam STO09_singles.sam
mv HTAdapter66_singles.sam STO10_singles.sam
mv HTAdapter67_singles.sam STO11_singles.sam
mv HTAdapter68_singles.sam STO12_singles.sam
mv HTAdapter69_singles.sam STO13_singles.sam
mv HTAdapter70_singles.sam STO14_singles.sam
mv HTAdapter71_singles.sam STO15_singles.sam
mv HTAdapter72_singles.sam STO17_singles.sam
mv HTAdapter73_singles.sam STO18_singles.sam
mv HTAdapter74_singles.sam STO19_singles.sam
mv HTAdapter75_singles.sam STO20_singles.sam
mv HTAdapter76_singles.sam STO21_singles.sam
mv HTAdapter77_singles.sam SUM01_singles.sam
mv HTAdapter78_singles.sam SUM02_singles.sam
mv HTAdapter79_singles.sam SUM03_singles.sam
mv HTAdapter80_singles.sam SUM04_singles.sam
mv HTAdapter81_singles.sam SUM05_singles.sam
mv HTAdapter82_singles.sam SUM06_singles.sam
mv HTAdapter83_singles.sam SUM07_singles.sam
mv HTAdapter84_singles.sam SUM08_singles.sam
mv HTAdapter85_singles.sam SUM09_singles.sam
mv HTAdapter86_singles.sam SUM10_singles.sam
mv HTAdapter87_singles.sam SUM11_singles.sam
mv HTAdapter88_singles.sam SUM12_singles.sam
mv HTAdapter89_singles.sam SUM13_singles.sam
mv HTAdapter90_singles.sam SUM15_singles.sam
mv HTAdapter91_singles.sam SUM16_singles.sam
mv HTAdapter92_singles.sam SUM17_singles.sam
mv HTAdapter93_singles.sam SUM18_singles.sam
mv HTAdapter94_singles.sam SUM19_singles.sam
mv HTAdapter95_singles.sam SUM20_singles.sam
mv HTAdapter96_singles.sam SUM21_singles.sam

# create list with each sam file
 ls | grep "_singles.sam" |cut -d "_" -f 1 | sort | uniq > readGroupList
 # while loop through the names in readGroupList
 while read i
 do
 # add read groups
 java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="$i"_singles.sam O="$i"_IDed.sam RGID="$i"_HTAdapter          RGLB=SeqCap2012 RGPL=illumina RGPU="$i" RGSM="$i"
 done<readGroupList
 echo "RLK_report: readGrouping complete"

# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# loop through .sam files
while read i;
do
# convert .sam to .bam & sort. ex: HTAdapter84_CACCTTAC_L003_R1.sam –> HTAdapter84_CACCTTAC_L003_R1_sorted.bam
samtools view -@ 2 -bS "$i"_IDed.sam | samtools sort -@ 2 -o "$i"_sorted.bam
# index the sorted .bam
samtools index "$i"_sorted.bam
# tally mapped reads & calcuate the stats
samtools idxstats "$i"_sorted.bam > "$i"_counts.txt
samtools flagstat "$i"_sorted.bam > "$i"_stats.txt
samtools depth ${wd}${i} > ${wd}${i}.depth.txt |awk '{sum+=$3} END {print sum/NR}' > "$i"_avgDepth.txt
done<samList

# make results directory & move results
mkdir -p /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/stats
mkdir -p /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/counts
mkdir -p /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/avgDepth

# cp results to respective directories
cp *stats.txt /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/stats
cp *counts.txt /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/counts
cp *depth.txt /home/rlk0015/SeqCap/Results_CompleteExonList/MapToExons/avgDepth

